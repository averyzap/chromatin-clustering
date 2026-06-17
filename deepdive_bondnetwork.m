%% DRIVER: bond-NETWORK deep-dive at a single lambda_cross
%  (companion to deepdive_lambda1.m, but cluster = bond-connected component)
%
%  Looks closely at one candidate flexible lambda_cross to see whether a
%  bond-connected group of >=3 beads PERSISTS while the individual edges
%  within it SWAP. That co-occurrence is the advisor's flexible regime.
%
%  Diagnostics:
%    1. Group size over time (size of largest bond-connected component)
%    2. Paired barcode: group membership (top) vs active edges (bottom)
%    3. "Edge swaps within a surviving group" count -- the direct flexible
%       signal: number of edge changes that happen WHILE the group membership
%       stays fixed
%    4. Replicate variability

clear; clc; close all;

%% ---------------- CONFIG ----------------
lambda_cross = 0.1;          % candidate flexible value; change after the sweep
sigma        = 60;
k_cross      = 0.1;

steps        = 300000;
warmup_steps = 30000;
n_replicates = 5;

min_size      = 3;
fine_stride_s = 0.02;

%% ---------------- RUN REPLICATES ----------------
reps = cell(1, n_replicates);
traces = cell(1, n_replicates);
for rep = 1:n_replicates
    fprintf('\n========== Replicate %d/%d (bond-network) ==========\n', rep, n_replicates);
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 600 + rep, 'save_traj', false, ...
                  'bond_stride', max(1, round(fine_stride_s / 1e-3)));
    reps{rep}   = run_hult_2d(lambda_cross, steps, warmup_steps, opts);
    traces{rep} = analyze_bondnetwork(reps{rep}, min_size);
end

%% ---------------- (1) GROUP SIZE TRACES ----------------
figure('Name', 'Bond-connected group size traces', 'Position', [50 50 1400 600])
for rep = 1:n_replicates
    subplot(n_replicates, 1, rep)
    tr = traces{rep};
    plot(tr.times, tr.group_size, 'b-', 'LineWidth', 0.8)
    ylim([0, max(2, max(tr.group_size)+0.5)])
    ylabel(sprintf('rep %d', rep))
    if rep == 1
        title(sprintf('Largest bond-connected group size over time (\\lambda_{cross} = %.3g)', lambda_cross))
    end
    if rep == n_replicates, xlabel('Time (s)'); end
    grid on
end

%% ---------------- (2) PAIRED BARCODE (replicate 1) ----------------
figure('Name', 'Bond-network paired barcode (rep 1)', 'Position', [100 100 1400 700])
tr = traces{1};
N  = tr.N;

subplot(2,1,1)
imagesc(tr.times, 1:N, double(tr.group_mat))
colormap(gca, [1 1 1; 0 0.4 0.8])
set(gca, 'YDir', 'reverse')
ylabel('bead #')
title(sprintf('Group membership (bond-connected \\geq%d)   \\lambda_{cross} = %.3g', min_size, lambda_cross))

subplot(2,1,2)
pair_labels = {}; pair_active = []; pidx = 0;
for i = 1:N
    for j = i+1:N
        if abs(i-j) > 1
            pidx = pidx + 1;
            pair_labels{pidx} = sprintf('%d-%d', i, j); %#ok<SAGROW>
            pair_active(pidx,:) = squeeze(tr.bond_mat(i,j,:)).'; %#ok<SAGROW>
        end
    end
end
imagesc(tr.times, 1:pidx, pair_active)
colormap(gca, [1 1 1; 0.8 0.2 0.2])
set(gca, 'YDir', 'reverse', 'YTick', 1:pidx, 'YTickLabel', pair_labels)
ylabel('bead pair'); xlabel('Time (s)')
title('Active edges (flexible = stable band above, swapping below)')

%% ---------------- (3) EDGE SWAPS WITHIN A SURVIVING GROUP ----------------
% The direct flexible-regime signal: count edge changes that occur while the
% group membership is unchanged. Many such swaps = flexible; few = rigid or
% amorphic.
swaps_within = zeros(1, n_replicates);
group_survival_time = zeros(1, n_replicates);
for rep = 1:n_replicates
    tr = traces{rep};
    swaps_within(rep) = tr.edge_swaps_within_group;
    group_survival_time(rep) = tr.median_group_survival;
end

figure('Name', 'Flexible signal: edge swaps within surviving groups', 'Position', [150 150 1000 450])
subplot(1,2,1)
bar(swaps_within, 'FaceColor', [0.2 0.5 0.8])
hold on; yline(mean(swaps_within), '--k', 'mean')
xlabel('Replicate'); ylabel('Edge swaps while group membership fixed')
title('Edge swaps within surviving groups')
grid on

subplot(1,2,2)
bar(group_survival_time, 'FaceColor', [0.2 0.4 0.8])
hold on; yline(median(group_survival_time), '--k', 'median')
xlabel('Replicate'); ylabel('Median group survival time (s)')
title('How long a bond-connected group survives')
grid on

%% ---------------- SUMMARY ----------------
fprintf('\n========== BOND-NETWORK DEEP-DIVE SUMMARY: lambda_cross = %.3g ==========\n', lambda_cross);
fprintf('%-8s %-18s %-22s %-22s\n', 'rep', 'mean group size', 'median group survival(s)', 'edge swaps within group');
fprintf('%s\n', repmat('-', 1, 75));
for rep = 1:n_replicates
    tr = traces{rep};
    fprintf('%-8d %-18.2f %-22.3f %-22d\n', rep, mean(tr.group_size), ...
        tr.median_group_survival, tr.edge_swaps_within_group);
end
fprintf('%s\n', repmat('-', 1, 75));
fprintf('Interpretation:\n');
fprintf('  Many edge swaps WHILE group survives  -> FLEXIBLE\n');
fprintf('  Group survives but few/no edge swaps  -> RIGID\n');
fprintf('  Group rarely survives at all          -> AMORPHIC\n');

save('deepdive_bondnetwork.mat', 'reps', 'traces', 'swaps_within', ...
     'group_survival_time', 'lambda_cross', 'sigma', 'k_cross', ...
     'min_size', 'fine_stride_s', 'n_replicates');
fprintf('\nSaved to deepdive_bondnetwork.mat\n');


%% =========================================================================
function tr = analyze_bondnetwork(r, min_size)
tr = struct();
N = r.N; tr.N = N;
if ~isfield(r, 'B_fine') || isempty(r.B_fine)
    error('analyze_bondnetwork requires opts.bond_stride > 0');
end
Bf = r.B_fine; times = r.B_fine_times(:).';
nfine = size(Bf, 3); dt_fine = times(1);

group_size = zeros(1, nfine);
group_mat  = false(N, nfine);
sets       = cell(1, nfine);
for k = 1:nfine
    Bk = (Bf(:,:,k) > 0);
    [lbl_of, comp_sizes] = connected_components(Bk);
    big = find(comp_sizes >= min_size);
    if isempty(big)
        sets{k} = []; group_size(k) = 0;
    else
        [gs, which] = max(comp_sizes(big));
        target = big(which);
        members = find(lbl_of == target);
        sets{k} = members(:); group_size(k) = gs;
        group_mat(members, k) = true;
    end
end

tr.times = times; tr.group_size = group_size;
tr.group_mat = group_mat; tr.bond_mat = Bf;

% group survival run-lengths (Jaccard > 0.5)
run_lengths = []; run_start = 1;
for k = 2:nfine
    if isempty(sets{k}) || isempty(sets{k-1})
        if ~isempty(sets{k-1}), run_lengths(end+1) = (k-1)-run_start+1; end %#ok<AGROW>
        run_start = k; continue
    end
    jac = numel(intersect(sets{k},sets{k-1}))/max(numel(union(sets{k},sets{k-1})),1);
    if jac < 0.5
        run_lengths(end+1) = (k-1)-run_start+1; %#ok<AGROW>
        run_start = k;
    end
end
if ~isempty(sets{nfine}), run_lengths(end+1) = nfine-run_start+1; end
if isempty(run_lengths), tr.median_group_survival = 0;
else, tr.median_group_survival = median(run_lengths)*dt_fine; end

% edge swaps that occur WHILE group membership is unchanged
swaps = 0;
for k = 2:nfine
    if isempty(sets{k}) || isempty(sets{k-1}), continue; end
    jac = numel(intersect(sets{k},sets{k-1}))/max(numel(union(sets{k},sets{k-1})),1);
    if jac >= 0.5
        % membership essentially unchanged; count edge changes among these beads
        members = union(sets{k}, sets{k-1});
        b1 = Bf(members, members, k-1);
        b2 = Bf(members, members, k);
        swaps = swaps + nnz(triu(xor(b1,b2),1));
    end
end
tr.edge_swaps_within_group = swaps;
end

%% =========================================================================
function [label, comp_sizes] = connected_components(A)
% Explicit connected-components via BFS on a logical symmetric adjacency
% matrix A. Every node gets a label; isolated nodes are singletons.
n = size(A, 1);
label = zeros(n, 1);
next_label = 0;
for s = 1:n
    if label(s) ~= 0, continue; end
    next_label = next_label + 1;
    queue = s;
    label(s) = next_label;
    while ~isempty(queue)
        v = queue(1); queue(1) = [];
        nbrs = find(A(v, :));
        for w = nbrs
            if label(w) == 0
                label(w) = next_label;
                queue(end+1) = w; %#ok<AGROW>
            end
        end
    end
end
comp_sizes = accumarray(label, 1);
end