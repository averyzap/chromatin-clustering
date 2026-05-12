%% TETHERED POLYMER SIMULATION (1D) -- HULT HARD-CUTOFF KERNEL

clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N     = 6;
L     = 120;
steps = 200000;
dt    = 1e-3;

zeta = 2.5e-3;
kBT  = 4.1;

% chain WLC
Lp    = 50;
Nk    = 17;
R0    = Nk*(2*Lp);     % chain WLC singularity
alpha = 0.2176;        % chain WLC stiffness coefficient

% Hult eligibility barrier
r_barrier = 90;            % nm; pairs with |x_i - x_j| < r_barrier are eligible

% Excluded volume
cEV = 8.305e-5;
aEV = 3.268e-5;

% Crosslink kinetics
kon0    = 1.0;             % on-rate inside the barrier (uniform)
koff    = 0.5;             % constant off-rate
k_cross = 0.01;            % linear crosslink spring stiffness (simplified)

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONDITION ----------------
x        = linspace(0, L, N)';
x(1)     = 0;
x(N)     = L;
B        = zeros(N);     % bond matrix (symmetric)
bound    = zeros(N, 1);  % per-bead bond flag
bond_start     = zeros(N);   % step at which each bond formed
bond_lifetimes = [];         % completed-bond lifetimes (s)
time_bonded    = zeros(N);   % cumulative steps each pair has been bonded

store_interval = 5000;
num_snaps      = floor(steps / store_interval);
B_history       = cell(1, num_snaps);
traj            = zeros(N, steps);
largest_cluster = zeros(num_snaps, 1);
cluster_sizes   = zeros(num_snaps, N);
snap_idx = 1;

%% ---------------- MAIN LOOP ----------------
fprintf('Running Hult hard-cutoff simulation (%d beads, %d steps)\n', N, steps);
tic
for t = 1:steps

    %% --- CROSSLINK DYNAMICS ---

    % Step 1: unbinding (each existing bond breaks independently with prob koff*dt)
    for i = 1:N
        for j = i+1:N
            if B(i,j) == 1 && rand < koff * dt
                lifetime               = (t - bond_start(i,j)) * dt;
                bond_lifetimes(end+1)  = lifetime; %#ok<AGROW>
                B(i,j)          = 0;
                B(j,i)          = 0;
                bound(i)        = 0;
                bound(j)        = 0;
                bond_start(i,j) = 0;
            end
        end
    end

    % Step 2: binding -- closest-pair priority within the eligibility barrier.
    % We build the list of eligible unbound pairs (|i-j|>1, both unbound,
    % within r_barrier), sort by distance, and attempt binding in order.
    eligible_pairs = [];
    eligible_dist  = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1 && B(i,j) == 0 && bound(i) == 0 && bound(j) == 0
                d = abs(x(i) - x(j));
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
            % re-check bound state (an earlier pair this step may have claimed i or j)
            if bound(i) == 0 && bound(j) == 0
                if rand < kon0 * dt
                    B(i,j)          = 1;
                    B(j,i)          = 1;
                    bound(i)        = 1;
                    bound(j)        = 1;
                    bond_start(i,j) = t;
                end
            end
        end
    end

    %% --- ACCUMULATE TIME BONDED ---
    time_bonded = time_bonded + B;

    %% --- FORCES ---
    F = zeros(N, 1);

    % chain WLC
    r    = diff(x);
    rabs = abs(r);
    FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);
    F(2:N)   = F(2:N)   - sign(r).*FWLC;
    F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

    % crosslink force (simplified linear spring)
    dx_cross = x' - x;
    Fcross   = k_cross * B .* dx_cross;
    F        = F + sum(Fcross, 2);

    % excluded volume
    dx_mat = x - x';
    FEV    = cEV * dx_mat .* exp(-aEV * dx_mat.^2);
    F      = F + sum(FEV, 2);

    % integrate
    xi = noise * randn(N, 1);
    x  = x + (dt/zeta)*F + xi;
    x(1) = 0;
    x(N) = L;

    traj(:,t) = x;

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
r_thresh = 20;
for k = 1:num_snaps
    x_snap = traj(:, k*store_interval);
    dx     = abs(x_snap - x_snap');
    A      = (dx < r_thresh) & (dx > 0);
    G      = graph(A);
    bins   = conncomp(G);
    sizes  = zeros(1, max(bins));
    for c = 1:max(bins)
        sizes(c) = sum(bins == c);
    end
    cluster_sizes(k, 1:length(sizes)) = sizes;
    largest_cluster(k) = max(sizes);
end

% bond_fraction(i,j): fraction of total simulation time pair (i,j) was bonded
bond_fraction = (time_bonded / 2) / steps;

%% ---------------- PLOTS ----------------
snap_times = (1:num_snaps) * store_interval * dt;

% --- Cluster growth ---
figure
plot(snap_times, largest_cluster, 'b', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Largest cluster size')
title('Cluster growth (Hult hard-cutoff)')
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

% --- Bond network snapshots ---
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_history{k})
    axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Bond network snapshots')

% --- Final configuration ---
figure; hold on
x_final = traj(:, end);
B_final = B_history{end};
plot(x_final, zeros(size(x_final)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
for i = 1:N
    for j = i+1:N
        if B_final(i,j) == 1
            plot([x_final(i), x_final(j)], [0, 0], 'b-', 'LineWidth', 2)
        end
    end
end
xlabel('Position (nm)')
title('Final configuration')
ylim([-1 1])
grid on

%% ---------------- SUMMARY ----------------
fprintf('\n========== HULT HARD-CUTOFF SUMMARY ==========\n')
if ~isempty(bond_lifetimes)
    eligible = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1
                eligible(end+1) = bond_fraction(i,j); %#ok<AGROW>
            end
        end
    end
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
fprintf('==============================================\n')