%% TETHERED POLYMER SIMULATION (2D) -- HULT HARD-CUTOFF KERNEL
%  WITH:
%  - tau_cross timescale separation between crosslink and polymer dynamics

clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N     = 6;                   % bead count (matches 1D version)
R_nuc = 175;                 % disk radius, nm -- scaled to fit a 6-bead chain
                             % (5 segments x ~70 nm at rest spans a 350 nm chord)
steps = 200000;
dt    = 1e-3;

zeta = 2.5e-3;
kBT  = 4.1;

% chain WLC
Lp    = 50;
Nk    = 17;
R0    = Nk*(2*Lp);           % chain WLC singularity, ~1700 nm (per spring)
alpha = 0.2176;              % chain WLC stiffness coefficient

% Hult eligibility barrier
r_barrier = 90;              % nm

% Excluded volume
cEV = 8.305e-5;
aEV = 3.268e-5;

% Crosslink kinetics
kon0    = 1.0;
koff    = 0.5;
k_cross = 0.01;              % linear crosslink spring stiffness

% Timescale separation between crosslink dynamics and polymer dynamics.
% tau_cross < 1  -> crosslinks faster than polymer (quasi-static crosslink limit, Anna)
% tau_cross = 1  -> baseline (matches earlier 1D runs)
% tau_cross > 1  -> crosslinks slower than polymer
% Both rates are scaled identically so the equilibrium bound fraction
% kon0/(kon0+koff) is preserved.
tau_cross = 0.01;
kon_eff   = kon0 / tau_cross;
koff_eff  = koff / tau_cross;

% Sanity check: keep per-step probabilities well below 1 so the
% rate*dt approximation is valid.
p_max = max(kon_eff, koff_eff) * dt;
if p_max > 0.1
    warning('tau_cross=%.3g gives per-step P=%.3f; consider reducing dt.', ...
            tau_cross, p_max);
end

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONFIGURATION (SAW between pinned ends) ----------------
% Pin endpoints to the boundary on opposite sides of the disk.
x = zeros(N, 2);
x(1, :) = [-R_nuc, 0];
x(N, :) = [ R_nuc, 0];

% Target step length ~ chain end-to-end distance / (N-1), softened so the
% chain has slack to fluctuate. We aim for steps of ~70 nm, well below R0.
target_step = 70;
min_separation = 15;         % minimum separation between non-bonded beads at init
                             % (scaled down from 30 nm to match smaller domain)

max_attempts = 5000;
success = false;
for attempt = 1:max_attempts
    x_try = x;
    ok = true;
    for i = 2:N-1
        % bias each step toward the far endpoint, but with random angular noise
        remaining = N - i + 1;
        target = (x(N,:) - x_try(i-1,:)) / remaining;
        target_dir = target / max(norm(target), eps);
        angle = atan2(target_dir(2), target_dir(1)) + 0.8 * (rand - 0.5) * pi;
        step = target_step * [cos(angle), sin(angle)];
        candidate = x_try(i-1,:) + step;
        % reject if outside disk or too close to existing beads
        if norm(candidate) > R_nuc - 10 || ...
                any(vecnorm(x_try(1:i-1,:) - candidate, 2, 2) < min_separation)
            ok = false;
            break
        end
        x_try(i,:) = candidate;
    end
    if ok
        % final reachability check: last interior bead must be within
        % a chain-step of the pinned endpoint
        if norm(x_try(N-1,:) - x(N,:)) < 3 * target_step
            x = x_try;
            success = true;
            break
        end
    end
end
if ~success
    warning('SAW init failed after %d attempts; using straight-line fallback', max_attempts);
    for i = 2:N-1
        x(i,:) = x(1,:) + (i-1)/(N-1) * (x(N,:) - x(1,:));
    end
end

%% ---------------- BOND BOOKKEEPING ----------------
B              = zeros(N);
bound          = zeros(N, 1);
bond_start     = zeros(N);
bond_lifetimes = [];
time_bonded    = zeros(N);

store_interval = 5000;
num_snaps      = floor(steps / store_interval);
B_history       = cell(1, num_snaps);
traj            = zeros(N, 2, steps);
largest_cluster = zeros(num_snaps, 1);
cluster_sizes   = zeros(num_snaps, N);
snap_idx = 1;

%% ---------------- MAIN LOOP ----------------
fprintf('Running 2D Hult hard-cutoff simulation (%d beads, %d steps, tau_cross=%.3g)\n', ...
        N, steps, tau_cross);
tic
for t = 1:steps

    %% --- CROSSLINK DYNAMICS ---

    % unbinding pass
    for i = 1:N
        for j = i+1:N
            if B(i,j) == 1 && rand < koff_eff * dt
                lifetime              = (t - bond_start(i,j)) * dt;
                bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
                B(i,j)          = 0;
                B(j,i)          = 0;
                bound(i)        = 0;
                bound(j)        = 0;
                bond_start(i,j) = 0;
            end
        end
    end

    % binding pass: closest-pair priority within eligibility barrier
    eligible_pairs = [];
    eligible_dist  = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1 && B(i,j) == 0 && bound(i) == 0 && bound(j) == 0
                d = norm(x(i,:) - x(j,:));
                if d < r_barrier
                    eligible_pairs(end+1, :) = [i, j]; %#ok<AGROW>
                    eligible_dist(end+1)     = d;     %#ok<AGROW>
                end
            end
        end
    end

    if ~isempty(eligible_dist)
        [~, order] = sort(eligible_dist);
        for k = 1:length(order)
            ij = eligible_pairs(order(k), :);
            i = ij(1); j = ij(2);
            if bound(i) == 0 && bound(j) == 0
                if rand < kon_eff * dt
                    B(i,j)          = 1;
                    B(j,i)          = 1;
                    bound(i)        = 1;
                    bound(j)        = 1;
                    bond_start(i,j) = t;
                end
            end
        end
    end

    % accumulate time bonded
    time_bonded = time_bonded + B;

    %% --- FORCES ---
    F = zeros(N, 2);

    % chain WLC (between consecutive beads)
    for i = 1:N-1
        rij  = x(i+1,:) - x(i,:);
        rabs = norm(rij);
        if rabs > 0 && rabs < R0
            f_mag = alpha * (-1 + 1/(1 - rabs/R0)^2 + 4*rabs/R0);
            f_vec = f_mag * rij / rabs;
            F(i,:)   = F(i,:)   + f_vec;
            F(i+1,:) = F(i+1,:) - f_vec;
        end
    end

    % crosslink force (simplified linear spring) -- vectorized over both dims
    for d = 1:2
        dx_d   = x(:,d).' - x(:,d);
        Fcross = k_cross * B .* dx_d;
        F(:,d) = F(:,d) + sum(Fcross, 2);
    end

    % excluded volume (Gaussian repulsion, applied component-wise)
    for d = 1:2
        dx_d = x(:,d) - x(:,d).';
        % use full pairwise squared distance for the exponent
        if d == 1
            r2 = (x(:,1) - x(:,1).').^2 + (x(:,2) - x(:,2).').^2;
        end
        FEV    = cEV * dx_d .* exp(-aEV * r2);
        F(:,d) = F(:,d) + sum(FEV, 2);
    end

    % integrate
    xi = noise * randn(N, 2);
    x  = x + (dt/zeta)*F + xi;

    % wall enforcement: project any bead outside the disk back to the boundary
    r_norms = vecnorm(x, 2, 2);
    outside = r_norms > R_nuc;
    if any(outside)
        x(outside, :) = R_nuc * x(outside, :) ./ r_norms(outside);
    end

    % re-pin the two endpoints
    x(1, :) = [-R_nuc, 0];
    x(N, :) = [ R_nuc, 0];

    traj(:,:,t) = x;

    if mod(t, store_interval) == 0
        B_history{snap_idx} = B;
        snap_idx = snap_idx + 1;
    end
end

% capture bonds still active at end of run
for i = 1:N
    for j = i+1:N
        if B(i,j) == 1
            lifetime              = (steps - bond_start(i,j)) * dt;
            bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
        end
    end
end

toc

%% ---------------- PROXIMITY CLUSTER ANALYSIS ----------------
r_thresh = 25;   % near the 1D value; scaled to the smaller domain
for k = 1:num_snaps
    x_snap = traj(:, :, k*store_interval);
    d_mat  = sqrt((x_snap(:,1) - x_snap(:,1)').^2 + (x_snap(:,2) - x_snap(:,2)').^2);
    A      = (d_mat < r_thresh) & (d_mat > 0);
    G      = graph(A);
    bins   = conncomp(G);
    sizes  = zeros(1, max(bins));
    for c = 1:max(bins)
        sizes(c) = sum(bins == c);
    end
    cluster_sizes(k, 1:length(sizes)) = sizes;
    largest_cluster(k) = max(sizes);
end

bond_fraction = (time_bonded / 2) / steps;

%% ---------------- PLOTS ----------------
snap_times = (1:num_snaps) * store_interval * dt;

% --- Initial vs final chain configuration ---
figure
subplot(1,2,1); hold on
theta_disk = linspace(0, 2*pi, 200);
plot(R_nuc*cos(theta_disk), R_nuc*sin(theta_disk), 'k-')
plot(traj(:,1,1), traj(:,2,1), 'b-o', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
plot(x(1,1), x(1,2), 'gs', 'MarkerSize', 12, 'LineWidth', 2)
plot(x(N,1), x(N,2), 'rs', 'MarkerSize', 12, 'LineWidth', 2)
axis equal; grid on
xlabel('x (nm)'); ylabel('y (nm)')
title('Initial configuration')

subplot(1,2,2); hold on
plot(R_nuc*cos(theta_disk), R_nuc*sin(theta_disk), 'k-')
x_final = traj(:,:,end);
plot(x_final(:,1), x_final(:,2), 'b-o', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
% draw active bonds
B_final = B_history{end};
for i = 1:N
    for j = i+1:N
        if B_final(i,j) == 1
            plot([x_final(i,1), x_final(j,1)], [x_final(i,2), x_final(j,2)], ...
                 'r-', 'LineWidth', 1.5)
        end
    end
end
plot(x_final(1,1), x_final(1,2), 'gs', 'MarkerSize', 12, 'LineWidth', 2)
plot(x_final(N,1), x_final(N,2), 'rs', 'MarkerSize', 12, 'LineWidth', 2)
axis equal; grid on
xlabel('x (nm)'); ylabel('y (nm)')
title('Final configuration (red = active crosslinks)')

% --- Cluster growth ---
figure
plot(snap_times, largest_cluster, 'b', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Largest cluster size')
title(sprintf('Cluster growth (Hult hard-cutoff, 2D, \\tau_{cross}=%.3g)', tau_cross))
grid on

% --- Cluster size distribution ---
figure
histogram(cluster_sizes(cluster_sizes > 0), 'BinEdges', 0.5:1:N+0.5, 'FaceColor', 'b')
xlabel('Cluster size')
ylabel('Frequency')
title('Cluster size distribution')
grid on

% --- Bond lifetime distribution ---
figure
if ~isempty(bond_lifetimes)
    histogram(bond_lifetimes, 30, 'FaceColor', 'b')
    xline(mean(bond_lifetimes), '--k', 'LineWidth', 1.5)
    xlabel('Bond lifetime (s)')
    ylabel('Count')
    title(sprintf('Bond lifetimes  |  mean = %.3f s, median = %.3f s', ...
        mean(bond_lifetimes), median(bond_lifetimes)))
else
    title('No bonds formed')
end
grid on

% --- Bond fraction heatmap ---
figure
bf = bond_fraction;
mask = abs((1:N)' - (1:N)) <= 1;
bf(mask) = NaN;
imagesc(bf)
colorbar
colormap('hot')
axis square
xlabel('Bead j'); ylabel('Bead i')
title('Fraction of time bonded per bead pair')

%% ---------------- SUMMARY ----------------
fprintf('\n========== 2D HULT HARD-CUTOFF SUMMARY ==========\n')
if ~isempty(bond_lifetimes)
    eligible = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1
                eligible(end+1) = bond_fraction(i,j); %#ok<AGROW>
            end
        end
    end
    fprintf('  N beads, R_nuc     : %d, %.0f nm\n', N, R_nuc)
    fprintf('  tau_cross          : %.3g\n', tau_cross)
    fprintf('  kon_eff, koff_eff  : %.3g, %.3g\n', kon_eff, koff_eff)
    fprintf('  Total bonds formed : %d\n',     length(bond_lifetimes))
    fprintf('  Mean lifetime      : %.4f s\n', mean(bond_lifetimes))
    fprintf('  Median lifetime    : %.4f s\n', median(bond_lifetimes))
    fprintf('  Max lifetime       : %.4f s\n', max(bond_lifetimes))
    fprintf('  Std deviation      : %.4f s\n', std(bond_lifetimes))
    fprintf('  Mean bond fraction : %.4f\n',   mean(eligible))
    fprintf('  Max bond fraction  : %.4f\n',   max(eligible))
else
    fprintf('  No bonds recorded\n')
end
fprintf('=================================================\n')