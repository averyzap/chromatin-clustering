%% TETHERED POLYMER SIMULATION (1D)
% WITH:
% - Single binding site per bead
% - Runs Caitlin (Gaussian) and Ben (logistic) kon models in parallel
% - Matched kon scaling at r = 31 nm
% - Proximity-based clustering
% - Bond lifetime tracking
% - Replication analysis (n_trials runs per model, mean +/- std)

clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N       = 6;
L       = 120;
steps   = 200000;
dt      = 1e-3;
n_trials = 10;          % number of independent replicates per model

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
xlabel('Distance (nm)'); ylabel('k_{on}')
legend('Gaussian (Caitlin)', 'Ben logistic (scaled)', 'Match point')
title('Matched k_{on} Functions'); grid on

%% ---------------- SETUP ----------------
model_names   = {'Caitlin (Gaussian)', 'Ben (Logistic)'};
use_ben_flags = [false, true];
colors        = {'b', 'r'};

store_interval = 5000;
num_snaps      = floor(steps / store_interval);

% Per-trial storage: rows = trials, cols = snapshots or pairs
% Scalar metrics collected across trials
for m = 1:2
    trial_stats(m).mean_lifetime   = zeros(n_trials, 1);
    trial_stats(m).median_lifetime = zeros(n_trials, 1);
    trial_stats(m).n_bonds         = zeros(n_trials, 1);
    trial_stats(m).mean_bf         = zeros(n_trials, 1);  % mean bond fraction
    trial_stats(m).max_bf          = zeros(n_trials, 1);  % max bond fraction
    trial_stats(m).mean_cluster    = zeros(n_trials, 1);  % mean largest cluster
    trial_stats(m).bond_fraction   = zeros(N, N, n_trials);
    trial_stats(m).largest_cluster = zeros(num_snaps, n_trials);
    % save full results only from last trial for detailed plots
    trial_stats(m).last_B_history  = {};
    trial_stats(m).last_traj       = [];
    trial_stats(m).last_lifetimes  = [];
    trial_stats(m).last_cluster_sizes = [];
end

%% ---------------- REPLICATION LOOP ----------------
for m = 1:2
    fprintf('\n=== Model: %s ===\n', model_names{m});
    use_ben = use_ben_flags(m);

    for tr = 1:n_trials
        fprintf('  Trial %d / %d\n', tr, n_trials);

        % Initial condition (random seed per trial)
        x     = linspace(0, L, N)';
        x(1)  = 0; x(N) = L;
        B        = zeros(N);
        bound    = zeros(N, 1);
        bond_start     = zeros(N);
        bond_lifetimes = [];
        time_bonded    = zeros(N);
        snap_idx = 1;
        B_history = cell(1, num_snaps);
        traj      = zeros(N, steps);
        largest_cluster = zeros(num_snaps, 1);
        cluster_sizes   = zeros(num_snaps, N);

        tic
        for t = 1:steps

            %% --- CROSSLINK DYNAMICS ---
            for i = 1:N
                for j = i+1:N
                    if abs(i-j) > 1

                        dist = x(i) - x(j);

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
                                B(i,j) = 1; B(j,i) = 1;
                                bound(i) = 1; bound(j) = 1;
                                bond_start(i,j) = t;
                            end

                        % unbinding
                        elseif B(i,j) == 1
                            if rand < koff * dt
                                lifetime = (t - bond_start(i,j)) * dt;
                                bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
                                B(i,j) = 0; B(j,i) = 0;
                                bound(i) = 0; bound(j) = 0;
                                bond_start(i,j) = 0;
                            end
                        end

                    end
                end
            end

            %% --- ACCUMULATE TIME BONDED ---
            time_bonded = time_bonded + B;

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
            x(1) = 0; x(N) = L;

            traj(:,t) = x;

            if mod(t, store_interval) == 0
                B_history{snap_idx} = B;
                snap_idx = snap_idx + 1;
            end

        end

        % capture bonds still active at end
        for i = 1:N
            for j = i+1:N
                if B(i,j) == 1
                    bond_lifetimes(end+1) = (steps - bond_start(i,j)) * dt; %#ok<AGROW>
                end
            end
        end
        toc

        %% --- PROXIMITY CLUSTER ANALYSIS ---
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

        %% --- STORE TRIAL METRICS ---
        bf = (time_bonded / 2) / steps;

        % collect eligible (non-adjacent) pair bond fractions
        eligible = [];
        for i = 1:N
            for j = i+1:N
                if abs(i-j) > 1
                    eligible(end+1) = bf(i,j); %#ok<AGROW>
                end
            end
        end

        trial_stats(m).mean_lifetime(tr)   = mean(bond_lifetimes);
        trial_stats(m).median_lifetime(tr) = median(bond_lifetimes);
        trial_stats(m).n_bonds(tr)         = length(bond_lifetimes);
        trial_stats(m).mean_bf(tr)         = mean(eligible);
        trial_stats(m).max_bf(tr)          = max(eligible);
        trial_stats(m).mean_cluster(tr)    = mean(largest_cluster);
        trial_stats(m).bond_fraction(:,:,tr) = bf;
        trial_stats(m).largest_cluster(:,tr) = largest_cluster;

        % keep last trial for detailed plots
        if tr == n_trials
            trial_stats(m).last_B_history     = B_history;
            trial_stats(m).last_traj          = traj;
            trial_stats(m).last_lifetimes     = bond_lifetimes;
            trial_stats(m).last_cluster_sizes = cluster_sizes;
        end

    end % trial loop
end % model loop

%% ---------------- SUMMARY TABLE ----------------
fprintf('\n')
fprintf('%-26s  %8s  %8s  %8s  %8s  %8s  %8s\n', ...
    'Model', 'N bonds', 'Mean LT', 'Med LT', 'Mean BF', 'Max BF', 'Mean Clust')
fprintf('%s\n', repmat('-', 1, 80))
for m = 1:2
    s = trial_stats(m);
    fprintf('%-26s  %5.1f±%-3.1f  %5.2f±%-4.2f  %5.2f±%-4.2f  %5.4f±%-6.4f  %5.4f±%-6.4f  %5.2f±%-4.2f\n', ...
        model_names{m}, ...
        mean(s.n_bonds),         std(s.n_bonds), ...
        mean(s.mean_lifetime),   std(s.mean_lifetime), ...
        mean(s.median_lifetime), std(s.median_lifetime), ...
        mean(s.mean_bf),         std(s.mean_bf), ...
        mean(s.max_bf),          std(s.max_bf), ...
        mean(s.mean_cluster),    std(s.mean_cluster))
end
fprintf('%s\n', repmat('-', 1, 80))

%% ---------------- REPLICATION PLOTS ----------------
snap_times = (1:num_snaps) * store_interval * dt;

%% --- Mean Bond Fraction: bar chart with error bars ---
figure
metrics     = {'mean_bf', 'max_bf'};
metric_lbls = {'Mean bond fraction', 'Max bond fraction'};
for p = 1:2
    subplot(1, 2, p)
    vals = [mean(trial_stats(1).(metrics{p})), mean(trial_stats(2).(metrics{p}))];
    errs = [std(trial_stats(1).(metrics{p})),  std(trial_stats(2).(metrics{p}))];
    b = bar(vals, 0.5, 'FaceColor', 'flat');
    b.CData = [0 0 1; 1 0 0];
    hold on
    errorbar(1:2, vals, errs, 'k.', 'LineWidth', 1.5, 'CapSize', 10)
    set(gca, 'XTickLabel', {'Caitlin', 'Ben'})
    ylabel(metric_lbls{p})
    title([metric_lbls{p}, sprintf('\n(mean ± std, n=%d)', n_trials)])
    grid on
end
sgtitle('Bond Fraction Comparison')

%% --- Mean Lifetime: bar chart with error bars ---
figure
metrics     = {'mean_lifetime', 'median_lifetime'};
metric_lbls = {'Mean lifetime (s)', 'Median lifetime (s)'};
for p = 1:2
    subplot(1, 2, p)
    vals = [mean(trial_stats(1).(metrics{p})), mean(trial_stats(2).(metrics{p}))];
    errs = [std(trial_stats(1).(metrics{p})),  std(trial_stats(2).(metrics{p}))];
    b = bar(vals, 0.5, 'FaceColor', 'flat');
    b.CData = [0 0 1; 1 0 0];
    hold on
    errorbar(1:2, vals, errs, 'k.', 'LineWidth', 1.5, 'CapSize', 10)
    set(gca, 'XTickLabel', {'Caitlin', 'Ben'})
    ylabel(metric_lbls{p})
    title([metric_lbls{p}, sprintf('\n(mean ± std, n=%d)', n_trials)])
    grid on
end
sgtitle('Bond Lifetime Comparison')

%% --- Mean Cluster Size over time with shaded std ---
figure; hold on
for m = 1:2
    lc  = trial_stats(m).largest_cluster;       % num_snaps x n_trials
    mu  = mean(lc, 2);
    sd  = std(lc, 0, 2);
    fill([snap_times, fliplr(snap_times)], ...
         [mu'+sd', fliplr(mu'-sd')], ...
         colors{m}, 'FaceAlpha', 0.2, 'EdgeColor', 'none')
    plot(snap_times, mu, colors{m}, 'LineWidth', 2)
end
xlabel('Time (s)'); ylabel('Largest cluster size')
title(sprintf('Cluster Growth (mean ± std, n=%d)', n_trials))
legend(model_names); grid on

%% --- Bond Fraction per trial (strip plot) ---
figure; hold on
jitter = 0.08;
for m = 1:2
    xs = m + jitter*(rand(n_trials,1) - 0.5);
    scatter(xs, trial_stats(m).mean_bf, 40, colors{m}, 'filled', 'MarkerFaceAlpha', 0.6)
    errorbar(m, mean(trial_stats(m).mean_bf), std(trial_stats(m).mean_bf), ...
        'k', 'LineWidth', 2, 'CapSize', 12)
end
xlim([0.5 2.5]); set(gca, 'XTick', [1 2], 'XTickLabel', {'Caitlin', 'Ben'})
ylabel('Mean bond fraction')
title(sprintf('Mean Bond Fraction per Trial (n=%d)', n_trials))
grid on

%% --- Mean Bond Fraction heatmap (averaged over trials) ---
figure
for m = 1:2
    subplot(1, 2, m)
    bf_mean = mean(trial_stats(m).bond_fraction, 3);
    mask = abs((1:N)' - (1:N)) <= 1;
    bf_mean(mask) = NaN;
    clim_max = max(max(mean(trial_stats(1).bond_fraction, 3), [], 'all'), ...
                   max(mean(trial_stats(2).bond_fraction, 3), [], 'all'));
    imagesc(bf_mean, [0, clim_max])
    colorbar; colormap('hot'); axis square
    xlabel('Bead j'); ylabel('Bead i')
    title(sprintf('%s\n(averaged over %d trials)', model_names{m}, n_trials))
end
sgtitle('Mean Fraction of Time Bonded per Bead Pair')

%% ---------------- SINGLE-TRIAL DETAIL PLOTS (last trial) ----------------

%% --- Bond Lifetime Distributions (last trial) ---
figure
for m = 1:2
    lt = trial_stats(m).last_lifetimes;
    subplot(1, 2, m)
    if ~isempty(lt)
        histogram(lt, 30, 'FaceColor', colors{m})
        xline(mean(lt), '--k', 'LineWidth', 1.5)
        xlabel('Bond lifetime (s)'); ylabel('Count')
        title(sprintf('%s\nMean=%.3fs, Median=%.3fs (trial %d)', ...
            model_names{m}, mean(lt), median(lt), n_trials))
    end
    grid on
end
sgtitle('Bond Lifetime Distributions (last trial)')

%% --- Bond Lifetime CDF (last trial) ---
figure; hold on
for m = 1:2
    lt = sort(trial_stats(m).last_lifetimes);
    if ~isempty(lt)
        plot(lt, (1:length(lt))/length(lt), colors{m}, 'LineWidth', 2)
    end
end
xlabel('Bond lifetime (s)'); ylabel('Cumulative probability')
title(sprintf('Bond Lifetime CDF (trial %d)', n_trials))
legend(model_names); grid on

%% --- Bond Matrices (last trial, Ben only) ---
B_hist = trial_stats(2).last_B_history;
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_hist{k}); axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle(['Bond Network — ', model_names{2}, ' (last trial)'])

%% --- Final Configurations (last trial) ---
figure
for m = 1:2
    subplot(1, 2, m); hold on
    x_final = trial_stats(m).last_traj(:, end);
    B_final = trial_stats(m).last_B_history{end};
    plot(x_final, zeros(size(x_final)), 'ko', 'MarkerFaceColor', 'k')
    for i = 1:N
        for j = i+1:N
            if B_final(i,j) == 1
                plot([x_final(i), x_final(j)], [0, 0], [colors{m}, '-'], 'LineWidth', 2)
            end
        end
    end
    title(['Final Config — ', model_names{m}, sprintf(' (trial %d)', n_trials)])
    ylim([-1 1]); grid on
end
sgtitle('Final Configurations (last trial)')