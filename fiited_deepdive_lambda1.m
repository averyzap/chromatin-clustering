%% DRIVER: deep-dive at lambda_cross = 1 (the flexible regime)
%
%  Goal: characterize the time-resolved dynamics at the parameter value
%  where every sweep observable hit its minimum and where the autocorrelation
%  crossover lives. A sweep gives you one number per lambda_cross. This
%  script asks what's *happening* at that one value.
%
%  Diagnostics (in order of increasing depth):
%    1. Cluster size time series           - is it hopping or flat?
%    2. Cluster-identity barcode plot      - which beads are in the
%                                            largest cluster, over time
%    3. Bond network barcode plot          - which pairs are bonded over
%                                            time; do the same pairs
%                                            keep re-bonding?
%    4. Inter-change time distribution     - distribution of the time
%                                            between cluster-identity
%                                            changes; exponential? bimodal?
%    5. Replicate variability              - 5 seeds; how much does the
%                                            qualitative picture vary
%                                            between independent runs?
%
%  Held fixed:
%    - lambda_cross = 1
%    - sigma = 60, k_cross = 0.1, N = 6, R_nuc = 175

clear; clc; close all;

%% ---------------- CONFIG ----------------
lambda_cross = 1;
sigma        = 60;
k_cross      = 0.1;

steps        = 300000;     % ~300 s — well longer than tau_p ~ 1s and tau_bonds at lambda=1
warmup_steps = 30000;

n_replicates = 5;

% Take much denser cluster snapshots than the default 5000. We want to
% resolve cluster-size hops well below the polymer-relaxation timescale so
% the persistence-time estimate is not quantized by the sampling grid.
% At dt = 1e-3, stride = 10 gives 0.01 s resolution (100x below tau_p ~ 1 s).
fine_snap_stride = 10;

%% ---------------- RUN REPLICATES ----------------
replicates = cell(1, n_replicates);
for rep = 1:n_replicates
    fprintf('\n========== Replicate %d/%d ==========\n', rep, n_replicates);
    opts = struct('sigma',     sigma, ...
                  'k_cross',   k_cross, ...
                  'seed',      300 + rep, ...
                  'save_traj', true, ...
                  'verbose',   true);
    replicates{rep} = run_hult_2d(lambda_cross, steps, warmup_steps, opts);
end

%% ---------------- FINE-GRAINED CLUSTER ID TRACE PER REPLICATE ----------------
%
% run_hult_2d's built-in cluster analysis runs every store_interval (5000)
% steps. We want denser. Walk the saved trajectory at fine_snap_stride and
% recompute the largest-cluster membership at each.

r_thresh = 25;
fine_traces = cell(1, n_replicates);

for rep = 1:n_replicates
    r = replicates{rep};
    nsteps = size(r.traj, 3);
    snap_indices = fine_snap_stride : fine_snap_stride : nsteps;
    n_fine = length(snap_indices);

    trace = struct();
    trace.times             = snap_indices * r.dt;
    trace.largest_size      = zeros(1, n_fine);
    trace.largest_members   = false(r.N, n_fine);  % rows = beads, cols = time
    trace.bonds_active      = false(r.N, r.N, n_fine);

    for k = 1:n_fine
        t_idx  = snap_indices(k);
        x_snap = r.traj(:, :, t_idx);
        d_mat  = sqrt((x_snap(:,1) - x_snap(:,1)').^2 + (x_snap(:,2) - x_snap(:,2)').^2);
        A      = (d_mat < r_thresh) & (d_mat > 0);
        G      = graph(A);
        bins   = conncomp(G);
        sizes  = accumarray(bins(:), 1);
        [biggest_size, biggest_label] = max(sizes);
        trace.largest_size(k)       = biggest_size;
        trace.largest_members(:, k) = (bins(:) == biggest_label);

        % Bond network at this snapshot: recover from B_history (nearest
        % stored snapshot in time). The original B_history is on the coarse
        % stride; for the bond barcode we'll interpolate by nearest neighbor.
    end

    % Reconstruct bond-network barcode at the SAME fine times by snapping
    % each fine time to the nearest coarse B_history snapshot. This is OK
    % because run_hult_2d already saved B every store_interval steps.
    coarse_store = round(r.snap_times(1) / r.dt);  % steps per coarse snapshot
    for k = 1:n_fine
        coarse_idx = max(1, round(snap_indices(k) / coarse_store));
        coarse_idx = min(coarse_idx, length(r.B_history));
        if ~isempty(r.B_history{coarse_idx})
            trace.bonds_active(:,:,k) = (r.B_history{coarse_idx} > 0);
        end
    end

    fine_traces{rep} = trace;
end

%% ---------------- DETECT CLUSTER-IDENTITY CHANGES ----------------
% At each fine snapshot, compare largest-cluster bead-set to the previous
% snapshot. Mark a "change event" when Jaccard < 0.5 (membership flipped
% by more than half).

change_events     = cell(1, n_replicates);
inter_change_dt   = cell(1, n_replicates);
persistence_each  = zeros(1, n_replicates);

for rep = 1:n_replicates
    trace = fine_traces{rep};
    n_fine = length(trace.times);
    is_change = false(1, n_fine);
    for k = 2:n_fine
        prev_set = find(trace.largest_members(:, k-1));
        curr_set = find(trace.largest_members(:, k));
        if isempty(prev_set) && isempty(curr_set)
            continue
        end
        inter = intersect(prev_set, curr_set);
        uni   = union(prev_set, curr_set);
        if isempty(uni)
            continue
        end
        jac = length(inter) / length(uni);
        if jac < 0.5
            is_change(k) = true;
        end
    end
    change_events{rep} = find(is_change);
    if numel(change_events{rep}) >= 2
        inter_change_dt{rep} = diff(trace.times(change_events{rep}));
        persistence_each(rep) = median(inter_change_dt{rep});
    else
        inter_change_dt{rep} = [];
        persistence_each(rep) = NaN;
    end
end

%% ---------------- PLOTS ----------------

% --- (1) Cluster size time series, all replicates ---
figure('Name', 'Cluster size traces at lambda_cross = 1', 'Position', [50 50 1400 600])
for rep = 1:n_replicates
    subplot(n_replicates, 1, rep)
    trace = fine_traces{rep};
    plot(trace.times, trace.largest_size, 'b-', 'LineWidth', 0.8)
    hold on
    if ~isempty(change_events{rep})
        for k = change_events{rep}
            xline(trace.times(k), 'r-', 'Alpha', 0.3);
        end
    end
    ylabel(sprintf('rep %d', rep))
    ylim([0.5, max(trace.largest_size) + 0.5])
    if rep == 1
        title(sprintf('Largest cluster size over time (\\lambda_{cross} = %.3g, %d replicates)', ...
                     lambda_cross, n_replicates))
    end
    if rep == n_replicates
        xlabel('Time (s)')
    end
    grid on
end

% --- (2) Cluster-identity barcode plots (one panel per replicate) ---
figure('Name', 'Cluster identity barcodes', 'Position', [100 100 1400 800])
for rep = 1:n_replicates
    subplot(n_replicates, 1, rep)
    trace = fine_traces{rep};
    imagesc(trace.times, 1:size(trace.largest_members, 1), double(trace.largest_members))
    colormap(gca, [1 1 1; 0 0.4 0.8])   % white = not in cluster, blue = in cluster
    set(gca, 'YDir', 'reverse')
    ylabel('bead #')
    if rep == 1
        title('Which beads are in the largest cluster, over time')
    end
    if rep == n_replicates
        xlabel('Time (s)')
    end
end

% --- (3) Bond network barcode, one representative replicate ---
figure('Name', 'Bond network over time (replicate 1)', 'Position', [150 150 1400 500])
trace = fine_traces{1};
N = size(trace.bonds_active, 1);
pair_idx = 0;
pair_labels = {};
pair_active = [];
for i = 1:N
    for j = i+1:N
        if abs(i-j) > 1   % only non-adjacent pairs are bondable
            pair_idx = pair_idx + 1;
            pair_labels{pair_idx} = sprintf('%d-%d', i, j);
            pair_active(pair_idx, :) = squeeze(trace.bonds_active(i, j, :)).';
        end
    end
end
imagesc(trace.times, 1:pair_idx, pair_active)
colormap([1 1 1; 0.8 0.2 0.2])
set(gca, 'YDir', 'reverse', 'YTick', 1:pair_idx, 'YTickLabel', pair_labels)
ylabel('bead pair')
xlabel('Time (s)')
title(sprintf('Bond network over time (replicate 1, \\lambda_{cross} = %.3g)', lambda_cross))

% --- (4) Inter-change time distribution, pooled across replicates ---
all_dt = [];
for rep = 1:n_replicates
    all_dt = [all_dt, inter_change_dt{rep}]; %#ok<AGROW>
end

figure('Name', 'Inter-change time distribution', 'Position', [200 200 900 500])
subplot(1,2,1)
if ~isempty(all_dt)
    histogram(all_dt, 20, 'FaceColor', [0.2 0.4 0.8])
    xlabel('Time between cluster-identity changes (s)')
    ylabel('Count')
    title(sprintf('Pooled inter-change time distribution (n=%d events)', length(all_dt)))
end
grid on

subplot(1,2,2)
tau_fit = NaN;
if ~isempty(all_dt) && length(all_dt) > 5
    dt_sorted = sort(all_dt, 'descend');
    surv      = (1:length(dt_sorted)) / length(dt_sorted);
    semilogy(dt_sorted, surv, 'b-', 'LineWidth', 1.5)
    hold on

    % Exponential fit: for a pure exponential, P(T>t) = exp(-t/tau), so
    % log(survival) is linear in t with slope -1/tau. Fit over the bulk
    % (survival between 0.01 and 0.9) to avoid the noisy extreme tail.
    t_asc    = sort(all_dt, 'ascend');
    surv_asc = 1 - (1:length(t_asc))'/length(t_asc);
    keep     = surv_asc > 0.01 & surv_asc < 0.9;
    if nnz(keep) > 5
        p = polyfit(t_asc(keep), log(surv_asc(keep)), 1);
        tau_fit = -1 / p(1);
        % overlay the fit
        tt = linspace(0, max(all_dt), 200);
        semilogy(tt, exp(p(2)) * exp(-tt / tau_fit), 'r--', 'LineWidth', 1.5)
        legend('data', sprintf('exp fit: \\tau = %.2f s', tau_fit), 'Location', 'best')
    end

    xlabel('Time between changes (s)')
    ylabel('Survival probability')
    title('Inter-change time survival (line straight = exponential)')
end
grid on

% --- (5) Replicate variability summary ---
figure('Name', 'Replicate variability', 'Position', [250 250 1000 400])
subplot(1,2,1)
mean_cluster_per_rep = arrayfun(@(rep) mean(fine_traces{rep}.largest_size), 1:n_replicates);
bar(mean_cluster_per_rep, 'FaceColor', [0.2 0.4 0.8])
hold on
yline(mean(mean_cluster_per_rep), '--k', 'mean across replicates')
xlabel('Replicate'); ylabel('Time-mean largest cluster')
title('Cluster size across replicates')
grid on

subplot(1,2,2)
bar(persistence_each, 'FaceColor', [0.6 0.2 0.6])
hold on
yline(median(persistence_each, 'omitnan'), '--k', 'median across replicates')
xlabel('Replicate'); ylabel('Median identity-persistence time (s)')
title('Persistence time across replicates')
grid on

%% ---------------- SUMMARY ----------------
fprintf('\n========== DEEP-DIVE SUMMARY: lambda_cross = %.3g ==========\n', lambda_cross);
fprintf('Replicates: %d   |   Run length per replicate: %.1f s   |   Fine snapshot stride: %.3f s\n', ...
        n_replicates, steps * replicates{1}.dt, fine_snap_stride * replicates{1}.dt);
fprintf('\n');
fprintf('%-8s %-18s %-22s %-18s\n', 'rep', '<largest cluster>', 'identity-changes count', 'persistence (s)');
fprintf('%s\n', repmat('-', 1, 75));
for rep = 1:n_replicates
    fprintf('%-8d %-18.2f %-22d %-18.3f\n', ...
        rep, mean(fine_traces{rep}.largest_size), ...
        length(change_events{rep}), persistence_each(rep));
end
fprintf('%s\n', repmat('-', 1, 75));
fprintf('Across all replicates: median persistence = %.3f s,  total change events = %d\n', ...
        median(persistence_each, 'omitnan'), length(all_dt));
if ~isempty(all_dt)
    fprintf('Inter-change time:   mean = %.3f s,  median = %.3f s,  max = %.3f s\n', ...
            mean(all_dt), median(all_dt), max(all_dt));
    if ~isnan(tau_fit)
        fprintf('Exponential fit to inter-change survival:  tau = %.3f s\n', tau_fit);
    end
end

save('deepdive_lambda1.mat', 'replicates', 'fine_traces', 'change_events', ...
     'inter_change_dt', 'persistence_each', 'all_dt', 'tau_fit', ...
     'lambda_cross', 'sigma', 'k_cross', 'n_replicates');
fprintf('\nSaved to deepdive_lambda1.mat\n');