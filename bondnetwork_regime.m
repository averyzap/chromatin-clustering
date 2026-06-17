%% DRIVER: bond-NETWORK regime analysis (companion to membership_bond_regime.m)
%
%  Same goal as membership_bond_regime.m -- identify the flexible regime --
%  but with cluster membership defined by the BOND NETWORK instead of spatial
%  proximity.
%
%    cluster = connected component of the bond graph B (beads linked through
%              active crosslinks), of size >= min_cluster_size
%
%  This sidesteps the proximity-flicker problem: a bond-connected group stays
%  defined as long as SOME bonds hold it together, even as individual bonds
%  swap. The advisor's regime definitions map cleanly onto two timescales:
%
%    RIGID    : same bond-connected group AND same edges       (both persist)
%    FLEXIBLE : same bond-connected group, edges keep swapping  (group persists,
%                                                                edges turn over)
%    AMORPHIC : bond-connected group itself dissolves           (group turns over)
%
%  NOTE on one-bond-per-bead: a bond-connected group of >=3 beads requires
%  >=2 simultaneous bonds sharing a bead (e.g. 3-4 and 4-5). Sustaining such a
%  group while edges swap requires a replacement bond to form before the old
%  one breaks. Whether that happens is exactly the empirical question here.

clear; clc; close all;

%% ---------------- CONFIG ----------------
lambda_cross_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5];
sigma             = 60;
k_cross           = 0.1;

steps_base        = 200000;
warmup_steps      = 20000;
steps_list        = round(steps_base * max(1, sqrt(lambda_cross_list)));

barcode_lambda    = 0.1;     % lambda_cross for the detailed paired-barcode figure

min_cluster_size  = 3;       % a bond-connected cluster must have >= this many beads
fine_stride_s     = 0.02;    % seconds between samples

n_sweep = length(lambda_cross_list);
results = cell(1, n_sweep);

%% ---------------- RUN SWEEP + ANALYZE ----------------
group_persist      = zeros(1, n_sweep);  % how long a >=3 bond-connected group survives (s)
edge_persist       = zeros(1, n_sweep);  % mean individual bond lifetime (s)
group_turnover     = zeros(1, n_sweep);  % fraction of group membership changing / s
edge_turnover      = zeros(1, n_sweep);  % fraction of edges changing / s
frac_time_grouped  = zeros(1, n_sweep);  % fraction of time a >=3 bond-group exists

barcode_trace = [];

fprintf('\n========== BOND-NETWORK REGIME SWEEP ==========\n');
fprintf('sigma = %.0f nm, k_cross = %.3g, min bond-group size = %d\n\n', ...
        sigma, k_cross, min_cluster_size);

for s = 1:n_sweep
    lc = lambda_cross_list(s);
    fprintf('--- %d/%d: lambda_cross = %.3g ---\n', s, n_sweep, lc);
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 500 + s, 'save_traj', false, ...
                  'bond_stride', max(1, round(fine_stride_s / 1e-3)));
    r = regime_analysis_hult_2d(lc, steps_list(s), warmup_steps, opts);

    trace = build_bondnetwork_trace(r, min_cluster_size);

    group_persist(s)     = trace.group_persist;
    edge_persist(s)      = trace.edge_persist;
    group_turnover(s)    = trace.group_turnover;
    edge_turnover(s)     = trace.edge_turnover;
    frac_time_grouped(s) = trace.frac_time_grouped;

    if abs(lc - barcode_lambda) < 1e-9
        barcode_trace = trace;
    end
    results{s} = r;
end

if isempty(barcode_trace)
    [~, idx] = min(abs(lambda_cross_list - barcode_lambda));
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 999, 'save_traj', false, ...
                  'bond_stride', max(1, round(fine_stride_s / 1e-3)));
    rb = regime_analysis_hult_2d(lambda_cross_list(idx), steps_list(idx), warmup_steps, opts);
    barcode_trace = build_bondnetwork_trace(rb, min_cluster_size);
    barcode_lambda = lambda_cross_list(idx);
end

%% ---------------- VISUAL 1: PERSISTENCE + TURNOVER ----------------
figure('Name', 'Bond-network: group vs edge', 'Position', [60 60 1300 550])

subplot(1,2,1)
loglog(lambda_cross_list, group_persist, 'bo-', 'LineWidth', 2, ...
       'MarkerFaceColor', 'b', 'DisplayName', 'Group persistence (same bond-connected beads)');
hold on
loglog(lambda_cross_list, edge_persist, 'rs-', 'LineWidth', 2, ...
       'MarkerFaceColor', 'r', 'DisplayName', 'Edge lifetime (individual bond, \propto \lambda_{cross})');
xlabel('\lambda_{cross}')
ylabel('Time (s)')
title({'Bond-connected group persistence vs edge lifetime', ...
       'Flexible: group persists while individual edges are short'})
legend('Location', 'best')
grid on

subplot(1,2,2)
semilogx(lambda_cross_list, group_turnover, 'bo-', 'LineWidth', 2, ...
         'MarkerFaceColor', 'b', 'DisplayName', 'Group membership turnover (/s)');
hold on
semilogx(lambda_cross_list, edge_turnover, 'rs-', 'LineWidth', 2, ...
         'MarkerFaceColor', 'r', 'DisplayName', 'Edge turnover (/s)');
xlabel('\lambda_{cross}')
ylabel('Turnover rate (fraction changing / s)')
title({'Turnover rates (same time windows)', ...
       'Flexible: edge turnover high while group turnover low'})
legend('Location', 'best')
grid on

%% ---------------- VISUAL 2: PAIRED BARCODE ----------------
figure('Name', sprintf('Bond-network barcode at lambda_cross = %.3g', barcode_lambda), ...
       'Position', [100 100 1400 700])

% Top: bond-connected group membership over time
subplot(2,1,1)
imagesc(barcode_trace.times, 1:barcode_trace.N, double(barcode_trace.group_mat))
colormap(gca, [1 1 1; 0 0.4 0.8])
set(gca, 'YDir', 'reverse')
ylabel('bead #')
title(sprintf(['Group: which beads are in the bond-connected (\\geq%d) cluster' ...
               '   (\\lambda_{cross} = %.3g)'], min_cluster_size, barcode_lambda))

% Bottom: which specific edges are active
subplot(2,1,2)
N = barcode_trace.N;
pair_labels = {};
pair_active = [];
pair_idx = 0;
for i = 1:N
    for j = i+1:N
        if abs(i-j) > 1
            pair_idx = pair_idx + 1;
            pair_labels{pair_idx} = sprintf('%d-%d', i, j); %#ok<SAGROW>
            pair_active(pair_idx, :) = squeeze(barcode_trace.bond_mat(i, j, :)).'; %#ok<SAGROW>
        end
    end
end
imagesc(barcode_trace.times, 1:pair_idx, pair_active)
colormap(gca, [1 1 1; 0.8 0.2 0.2])
set(gca, 'YDir', 'reverse', 'YTick', 1:pair_idx, 'YTickLabel', pair_labels)
ylabel('bead pair')
xlabel('Time (s)')
title('Edges: which specific bonds are active (flexible = stable group above, swapping edges below)')

%% ---------------- VISUAL 3: TURNOVER QUADRANT ----------------
figure('Name', 'Bond-network regime quadrant', 'Position', [150 150 800 650])
scatter(edge_turnover, group_turnover, 120, log10(lambda_cross_list), 'filled')
hold on
for s = 1:n_sweep
    text(edge_turnover(s)*1.05, group_turnover(s), ...
         sprintf('%.2g', lambda_cross_list(s)), 'FontSize', 9);
end
xlabel('Edge turnover (fraction of bonds changing / s)')
ylabel('Group turnover (fraction of bond-group membership changing / s)')
title({'Bond-network regime quadrant', ...
       'RIGID: low-low | FLEXIBLE: high edge turnover, low group turnover | AMORPHIC: high group turnover'})
cb = colorbar; cb.Label.String = 'log_{10}(\lambda_{cross})';
grid on

%% ---------------- SUMMARY ----------------
fprintf('\n========== BOND-NETWORK REGIME SUMMARY ==========\n');
fprintf('%-12s %-16s %-14s %-18s %-16s %-14s\n', ...
    'lambda', 'group_persist(s)', 'edge_persist(s)', 'group_turnover(/s)', 'edge_turnover(/s)', 'frac_grouped');
fprintf('%s\n', repmat('-', 1, 95));
for s = 1:n_sweep
    fprintf('%-12.3g %-16.3f %-14.3f %-18.3f %-16.3f %-14.3f\n', ...
        lambda_cross_list(s), group_persist(s), edge_persist(s), ...
        group_turnover(s), edge_turnover(s), frac_time_grouped(s));
end
fprintf('\nFlexible-regime signal = group_persist / edge_persist >> 1:\n');
ratio = group_persist ./ max(edge_persist, eps);
for s = 1:n_sweep
    fprintf('  lambda = %-8.3g  group/edge ratio = %.2f   (frac time grouped = %.2f)\n', ...
        lambda_cross_list(s), ratio(s), frac_time_grouped(s));
end

save('bondnetwork_regime.mat', 'lambda_cross_list', 'group_persist', ...
     'edge_persist', 'group_turnover', 'edge_turnover', 'frac_time_grouped', ...
     'barcode_trace', 'barcode_lambda', 'sigma', 'k_cross', ...
     'min_cluster_size', 'fine_stride_s');
fprintf('\nSaved to bondnetwork_regime.mat\n');


%% =========================================================================
function trace = build_bondnetwork_trace(r, min_size)
% Walk the fine bond storage. At each sample:
%   - build the bond graph from B_fine
%   - find the largest bond-connected component of size >= min_size (the group)
% Then compute group persistence, edge persistence, and turnover rates.

trace = struct();
N = r.N;
trace.N = N;

if ~isfield(r, 'B_fine') || isempty(r.B_fine)
    error('build_bondnetwork_trace requires opts.bond_stride > 0 (B_fine missing)');
end

Bf      = r.B_fine;               % N x N x nfine logical
times   = r.B_fine_times(:).';
nfine   = size(Bf, 3);
dt_fine = times(1);               % first timestamp = one stride in seconds

group_mat = false(N, nfine);
sets      = cell(1, nfine);

for k = 1:nfine
    Bk = (Bf(:, :, k) > 0);          % ensure logical adjacency
    [lbl_of, comp_sizes] = connected_components(Bk);
    big = find(comp_sizes >= min_size);
    if isempty(big)
        sets{k} = [];
    else
        [~, which] = max(comp_sizes(big));
        target = big(which);
        members = find(lbl_of == target);
        sets{k} = members(:);
        group_mat(members, k) = true;
    end
end

trace.times     = times;
trace.group_mat = group_mat;
trace.bond_mat  = Bf;

% --- group persistence: run-length of stable >=3 bond-connected set ---
run_lengths = [];
run_start = 1;
for k = 2:nfine
    if isempty(sets{k}) || isempty(sets{k-1})
        if ~isempty(sets{k-1})
            run_lengths(end+1) = (k-1) - run_start + 1; %#ok<AGROW>
        end
        run_start = k;
        continue;
    end
    inter = intersect(sets{k}, sets{k-1});
    uni   = union(sets{k}, sets{k-1});
    jac   = numel(inter) / max(numel(uni), 1);
    if jac < 0.5
        run_lengths(end+1) = (k-1) - run_start + 1; %#ok<AGROW>
        run_start = k;
    end
end
if ~isempty(sets{nfine})
    run_lengths(end+1) = nfine - run_start + 1;
end
if isempty(run_lengths)
    trace.group_persist = 0;
else
    trace.group_persist = median(run_lengths) * dt_fine;
end

% --- edge persistence: mean individual bond lifetime ---
if ~isempty(r.bond_lifetimes)
    trace.edge_persist = mean(r.bond_lifetimes);
else
    trace.edge_persist = 0;
end

% --- group membership turnover (Jaccard distance / s) ---
gch = [];
for k = 2:nfine
    if isempty(sets{k}) && isempty(sets{k-1}), continue; end
    inter = intersect(sets{k}, sets{k-1});
    uni   = union(sets{k}, sets{k-1});
    gch(end+1) = 1 - numel(inter)/max(numel(uni),1); %#ok<AGROW>
end
if isempty(gch), trace.group_turnover = 0; else, trace.group_turnover = mean(gch)/dt_fine; end

% --- edge turnover (fraction of bond set changing / s) ---
ech = [];
for k = 2:nfine
    b1 = Bf(:,:,k-1); b2 = Bf(:,:,k);
    n1 = nnz(triu(b1,1)); n2 = nnz(triu(b2,1));
    if n1==0 && n2==0, continue; end
    changed = nnz(triu(xor(b1,b2),1));
    ech(end+1) = changed / max(n1+n2,1); %#ok<AGROW>
end
if isempty(ech), trace.edge_turnover = 0; else, trace.edge_turnover = mean(ech)/dt_fine; end

% --- fraction of time a >=min_size bond-connected group exists ---
trace.frac_time_grouped = mean(cellfun(@(s) ~isempty(s), sets));
end

%% =========================================================================
function [label, comp_sizes] = connected_components(A)
% Explicit connected-components via BFS on a logical symmetric adjacency
% matrix A (no dependence on graph()/conncomp()). Every node gets a label;
% isolated nodes are their own singleton component.
n = size(A, 1);
label = zeros(n, 1);
next_label = 0;
for s = 1:n
    if label(s) ~= 0, continue; end
    next_label = next_label + 1;
    % BFS from s
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