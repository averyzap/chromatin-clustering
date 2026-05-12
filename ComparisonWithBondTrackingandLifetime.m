%% TETHERED POLYMER SIMULATION (1D)
% WITH:
% - Single binding site per bead
% - Runs Caitlin (Gaussian) and Ben (logistic) kon models in parallel
% - Matched kon scaling at r = 31 nm
% - Proximity-based clustering
% - Bond lifetime tracking

clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N = 6;
L = 120;
steps = 200000;
dt = 1e-3;

zeta = 2.5e-3;
kBT  = 4.1;

Lp = 50;
Nk = 17;
R0 = Nk*(2*Lp);
alpha = 0.2176;

% Excluded volume
cEV = 8.305e-5;
aEV = 3.268e-5;

% Crosslink parameters
kon0    = 1.0;
koff    = 0.5;
k_cross = 0.01;

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- BEN KERNEL MATCHING ----------------
r_match = 31;
r_half  = 21.5;
steep   = 5;

kG           = kon0 * exp(-(r_match^2)/(20^2));
a_match      = 2 / (1 + exp(steep*(r_match/r_half - 1)));
scale_factor = kG / (kon0 * a_match);

fprintf('Ben model scale factor = %f\n', scale_factor);

%% ---------------- VISUALIZE KON FUNCTIONS ----------------
r_vals  = linspace(0, 80, 500);
k_gauss = kon0 * exp(-(r_vals.^2)/(20^2));
k_ben   = zeros(size(r_vals));
for i = 1:length(r_vals)
    a        = 2 / (1 + exp(steep*(r_vals(i)/r_half - 1)));
    k_ben(i) = scale_factor * kon0 * a;
end

figure
plot(r_vals, k_gauss, 'b', 'LineWidth', 2); hold on
plot(r_vals, k_ben,   'r', 'LineWidth', 2)
xline(r_match, '--k', 'LineWidth', 1)
xlabel('Distance (nm)')
ylabel('k_{on}')
legend('Gaussian (Caitlin)', 'Ben logistic (scaled)', 'Match point')
title('Matched k_{on} Functions')
grid on

%% ---------------- RUN BOTH MODELS ----------------
model_names   = {'Caitlin (Gaussian)', 'Ben (Logistic)'};
use_ben_flags = [false, true];
colors        = {'b', 'r'};

store_interval = 5000;
num_snaps      = floor(steps / store_interval);

% Pre-allocate results struct
for m = 1:2
    results(m).B_history       = cell(1, num_snaps);
    results(m).traj            = zeros(N, steps);
    results(m).largest_cluster = zeros(num_snaps, 1);
    results(m).cluster_sizes   = zeros(num_snaps, N);
    results(m).bond_lifetimes  = [];
    results(m).time_bonded     = zeros(N);   % steps each pair spent bonded
end

for m = 1:2

    fprintf('\nRunning model: %s\n', model_names{m});
    use_ben = use_ben_flags(m);

    % Initial condition
    x        = linspace(0, L, N)';
    x(1)     = 0;
    x(N)     = L;
    B        = zeros(N);
    bound    = zeros(N, 1);
    bond_start     = zeros(N);   % step when each bond formed
    bond_lifetimes = [];
    time_bonded    = zeros(N);   % cumulative steps each pair is bonded
    snap_idx = 1;

    tic
    for t = 1:steps

        %% --- CROSSLINK DYNAMICS ---
        for i = 1:N
            for j = i+1:N
                if abs(i-j) > 1

                    dist = x(i) - x(j);

                    % compute kon
                    if use_ben
                        a      = 2 / (1 + exp(steep*(abs(dist)/r_half - 1)));
                        kon_ij = scale_factor * kon0 * a;
                    else
                        kon_ij = kon0 * exp(-(dist^2)/(20^2));
                    end
                    kon_ij = max(0, min(kon_ij, 1/dt));

                    % binding
                    if B(i,j) == 0 && bound(i) == 0 && bound(j) == 0
                        if rand < kon_ij * dt
                            B(i,j)          = 1;
                            B(j,i)          = 1;
                            bound(i)        = 1;
                            bound(j)        = 1;
                            bond_start(i,j) = t;
                        end

                    % unbinding
                    elseif B(i,j) == 1
                        if rand < koff * dt
                            lifetime             = (t - bond_start(i,j)) * dt;
                            bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
                            B(i,j)          = 0;
                            B(j,i)          = 0;
                            bound(i)        = 0;
                            bound(j)        = 0;
                            bond_start(i,j) = 0;
                        end
                    end

                end
            end
        end

        %% --- ACCUMULATE TIME BONDED ---
        time_bonded = time_bonded + B;   % B is symmetric; we'll halve later

        %% --- FORCES ---
        F    = zeros(N, 1);
        r    = diff(x);
        rabs = abs(r);
        FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);
        F(2:N)   = F(2:N)   - sign(r).*FWLC;
        F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

        dx     = x' - x;
        Fcross = k_cross * B .* dx;
        F      = F + sum(Fcross, 2);

        dx  = x - x';
        FEV = cEV * dx .* exp(-aEV * dx.^2);
        F   = F + sum(FEV, 2);

        xi = noise * randn(N, 1);
        x  = x + (dt/zeta)*F + xi;
        x(1) = 0;
        x(N) = L;

        results(m).traj(:,t) = x;

        if mod(t, store_interval) == 0
            results(m).B_history{snap_idx} = B;
            snap_idx = snap_idx + 1;
        end

    end

    % capture any bonds still active at end of run
    for i = 1:N
        for j = i+1:N
            if B(i,j) == 1
                lifetime             = (steps - bond_start(i,j)) * dt;
                bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
            end
        end
    end

    toc

    %% --- PROXIMITY CLUSTER ANALYSIS ---
    r_thresh = 20;
    for k = 1:num_snaps
        x_snap = results(m).traj(:, k*store_interval);
        dx     = abs(x_snap - x_snap');
        A      = (dx < r_thresh) & (dx > 0);
        G      = graph(A);
        bins   = conncomp(G);
        sizes  = zeros(1, max(bins));
        for c = 1:max(bins)
            sizes(c) = sum(bins == c);
        end
        results(m).cluster_sizes(k, 1:length(sizes)) = sizes;
        results(m).largest_cluster(k) = max(sizes);
    end

    results(m).bond_lifetimes = bond_lifetimes;

    % bond_fraction(i,j): fraction of total simulation time pair (i,j) was bonded
    % time_bonded is symmetric so divide by 2 to avoid double-counting,
    % then divide by total steps to get fraction in [0,1]
    results(m).bond_fraction = (time_bonded / 2) / steps;

end

%% ---------------- COMPARISON PLOTS ----------------
snap_times = (1:num_snaps) * store_interval * dt;

%% --- Cluster Growth ---
figure
for m = 1:2
    plot(snap_times, results(m).largest_cluster, colors{m}, 'LineWidth', 2); hold on
end
xlabel('Time (s)')
ylabel('Largest cluster size')
title('Cluster Growth Comparison')
legend(model_names)
grid on

%% --- Cluster Size Distribution ---
figure
for m = 1:2
    subplot(1, 2, m)
    sz = results(m).cluster_sizes;
    histogram(sz(sz > 0), 'BinEdges', 0.5:1:N+0.5, 'FaceColor', colors{m})
    xlabel('Cluster size')
    ylabel('Frequency')
    title(model_names{m})
    grid on
end
sgtitle('Cluster Size Distribution')

%% --- Bond Lifetime Distributions ---
figure
for m = 1:2
    lt = results(m).bond_lifetimes;
    subplot(1, 2, m)
    if ~isempty(lt)
        histogram(lt, 30, 'FaceColor', colors{m})
        xline(mean(lt), '--k', 'LineWidth', 1.5)
        xlabel('Bond lifetime (s)')
        ylabel('Count')
        title(sprintf('%s\nMean = %.3f s,  Median = %.3f s', ...
            model_names{m}, mean(lt), median(lt)))
    else
        title([model_names{m}, ' — no bonds formed'])
    end
    grid on
end
sgtitle('Bond Lifetime Distributions')

%% --- Bond Lifetime CDF (overlay) ---
figure; hold on
for m = 1:2
    lt = sort(results(m).bond_lifetimes);
    if ~isempty(lt)
        cdf = (1:length(lt)) / length(lt);
        plot(lt, cdf, colors{m}, 'LineWidth', 2)
    end
end
xlabel('Bond lifetime (s)')
ylabel('Cumulative probability')
title('Bond Lifetime CDF')
legend(model_names)
grid on

%% --- Bond Fraction Heatmaps ---
figure
for m = 1:2
    subplot(1, 2, m)
    bf = results(m).bond_fraction;
    % mask diagonal and nearest-neighbour bonds (not eligible for crosslinks)
    mask = abs((1:N)' - (1:N)) <= 1;
    bf(mask) = NaN;
    imagesc(bf, [0, max(results(1).bond_fraction(:), [], 'omitnan')])
    colorbar
    colormap('hot')
    axis square
    xlabel('Bead j'); ylabel('Bead i')
    title(model_names{m})
end
sgtitle('Fraction of Time Bonded per Bead Pair')


fprintf('\n========== BOND LIFETIME SUMMARY ==========\n')
for m = 1:2
    lt = results(m).bond_lifetimes;
    bf = results(m).bond_fraction;
    eligible = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1
                eligible(end+1) = bf(i,j); %#ok<AGROW>
            end
        end
    end
    if ~isempty(lt)
        fprintf('%s:\n',                        model_names{m})
        fprintf('  Total bonds formed : %d\n',     length(lt))
        fprintf('  Mean lifetime      : %.4f s\n', mean(lt))
        fprintf('  Median lifetime    : %.4f s\n', median(lt))
        fprintf('  Max lifetime       : %.4f s\n', max(lt))
        fprintf('  Std deviation      : %.4f s\n', std(lt))
        fprintf('  Mean bond fraction : %.4f\n',   mean(eligible))
        fprintf('  Max bond fraction  : %.4f\n',   max(eligible))
    else
        fprintf('%s: no bonds recorded\n', model_names{m})
    end
end
fprintf('============================================\n')

%% --- Bond Matrices ---
for m = 1:2
    figure
    for k = 1:num_snaps
        subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
        imagesc(results(m).B_history{k})
        axis square
        title(['t = ', num2str(k*store_interval)])
    end
    sgtitle(['Bond Network — ', model_names{m}])
end

%% --- Final Configurations ---
figure
for m = 1:2
    subplot(1, 2, m); hold on
    x_final = results(m).traj(:, end);
    B_final = results(m).B_history{end};
    plot(x_final, zeros(size(x_final)), 'ko', 'MarkerFaceColor', 'k')
    for i = 1:N
        for j = i+1:N
            if B_final(i,j) == 1
                plot([x_final(i), x_final(j)], [0, 0], [colors{m}, '-'], 'LineWidth', 2)
            end
        end
    end
    title(['Final Config — ', model_names{m}])
    ylim([-1 1]); grid on
end
sgtitle('Final Configurations')