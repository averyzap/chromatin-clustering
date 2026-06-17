%% DRIVER: membership-vs-bond regime analysis
%
%  Identifies the FLEXIBLE regime using the advisor's definition:
%
%    RIGID    : same beads stay clustered AND same bonds hold them
%               -> membership persistence LONG, bond persistence LONG
%    FLEXIBLE : same beads stay clustered, but the specific bonds among
%               them keep swapping
%               -> membership persistence LONG, bond persistence SHORT
%    AMORPHIC : the beads themselves don't stay clustered
%               -> membership persistence SHORT
%
%  The flexible regime is therefore a SEPARATION of two timescales:
%  spatial-cluster membership (long) vs bond-network within the cluster
%  (short). A single persistence metric cannot see this; we measure BOTH.
%
%  Cluster definition:
%    - MEMBERSHIP = proximity graph connected component of size >= 3
%      (any persistent group of >=3 beads, not necessarily the largest)
%    - BONDS      = the crosslink B matrix restricted to the clustered beads
%
%  Produces three visuals:
%    1. Two persistence curves vs lambda_cross (membership vs bond) -- the
%       flexible regime is where they SPLIT.
%    2. Paired membership/bond barcode at a candidate flexible lambda_cross.
%    3. Turnover scatter: bond turnover (x) vs membership turnover (y),
%       one point per lambda_cross, placing each into a regime quadrant.

clear; clc; close all;

%% ---------------- CONFIG ----------------
% Dense coverage between 0.01 and 1 (where the rigid->flexible transition
% is expected), plus a couple of points above 1 for the amorphic side.
lambda_cross_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5];
sigma             = 60;
k_cross           = 0.1;

steps_base        = 200000;
warmup_steps      = 20000;
steps_list        = round(steps_base * max(1, sqrt(lambda_cross_list)));

% Which lambda_cross to use for the detailed paired-barcode figure (visual 2).
% Default: the value you most suspect is flexible. Change after seeing visual 1.
barcode_lambda    = 0.03;

% Cluster / analysis parameters
r_thresh          = 25;       % nm, proximity threshold for "clustered"
min_cluster_size  = 3;        % a "cluster" must have >= this many beads
fine_stride_s     = 0.02;     % seconds between membership/bond samples
jaccard_thresh    = 0.5;      % membership "survives" if overlap > this

n_sweep = length(lambda_cross_list);
results = cell(1, n_sweep);

%% ---------------- RUN SWEEP + ANALYZE ----------------
membership_persist = zeros(1, n_sweep);   % how long a >=3 group stays together (s)
bond_persist       = zeros(1, n_sweep);   % how long an individual bond lasts (s)
membership_turnover= zeros(1, n_sweep);   % fraction of membership changing per second
bond_turnover      = zeros(1, n_sweep);   % fraction of bonds changing per second
frac_time_clustered= zeros(1, n_sweep);   % fraction of time a >=3 cluster exists at all

% store the trace for the barcode lambda
barcode_trace = [];

fprintf('\n========== MEMBERSHIP-vs-BOND REGIME SWEEP ==========\n');
fprintf('sigma = %.0f nm, k_cross = %.3g, r_thresh = %.0f nm, min size = %d\n\n', ...
        sigma, k_cross, r_thresh, min_cluster_size);

for s = 1:n_sweep
    lc = lambda_cross_list(s);
    fprintf('--- %d/%d: lambda_cross = %.3g ---\n', s, n_sweep, lc);
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 400 + s, 'save_traj', true, ...
                  'bond_stride', max(1, round(fine_stride_s / 1e-3)));
    r = regime_analysis_hult_2d(lc, steps_list(s), warmup_steps, opts);

    trace = build_membership_bond_trace(r, r_thresh, min_cluster_size, fine_stride_s);

    membership_persist(s)  = trace.membership_persist;
    bond_persist(s)        = trace.bond_persist;
    membership_turnover(s) = trace.membership_turnover;
    bond_turnover(s)       = trace.bond_turnover;
    frac_time_clustered(s) = trace.frac_time_clustered;

    % keep the full trace for the chosen barcode lambda
    if abs(lc - barcode_lambda) < 1e-9
        barcode_trace = trace;
    end

    % strip trajectory to save memory
    r.traj = [];
    results{s} = r;
end

% Fallback: if barcode_lambda wasn't in the list, use the closest run.
if isempty(barcode_trace)
    [~, idx] = min(abs(lambda_cross_list - barcode_lambda));
    fprintf('barcode_lambda %.3g not in sweep; using closest = %.3g\n', ...
            barcode_lambda, lambda_cross_list(idx));
    % rerun just that one to get its trace
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 999, 'save_traj', true, ...
                  'bond_stride', max(1, round(fine_stride_s / 1e-3)));
    rb = regime_analysis_hult_2d(lambda_cross_list(idx), steps_list(idx), warmup_steps, opts);
    barcode_trace = build_membership_bond_trace(rb, r_thresh, min_cluster_size, fine_stride_s);
    barcode_lambda = lambda_cross_list(idx);
end

%% ---------------- VISUAL 1: TWO PERSISTENCE CURVES + RATIO ----------------
figure('Name', 'Membership vs bond persistence', 'Position', [60 60 1300 550])

subplot(1,2,1)
loglog(lambda_cross_list, membership_persist, 'bo-', 'LineWidth', 2, ...
       'MarkerFaceColor', 'b', 'DisplayName', 'Membership persistence (same beads clustered)');
hold on
loglog(lambda_cross_list, bond_persist, 'rs-', 'LineWidth', 2, ...
       'MarkerFaceColor', 'r', 'DisplayName', 'Bond lifetime (global, \propto \lambda_{cross})');
xlabel('\lambda_{cross}')
ylabel('Time (s)')
title({'Membership persistence vs bond lifetime', ...
       'Flexible regime: membership stays high while bonds are short'})
legend('Location', 'best')
grid on

subplot(1,2,2)
% The real flexible-regime signal: membership turnover vs bond turnover,
% both measured over the SAME fine windows (independent of global lifetime).
semilogx(lambda_cross_list, membership_turnover, 'bo-', 'LineWidth', 2, ...
         'MarkerFaceColor', 'b', 'DisplayName', 'Membership turnover (/s)');
hold on
semilogx(lambda_cross_list, bond_turnover, 'rs-', 'LineWidth', 2, ...
         'MarkerFaceColor', 'r', 'DisplayName', 'Bond turnover (/s)');
xlabel('\lambda_{cross}')
ylabel('Turnover rate (fraction changing / s)')
title({'Turnover rates (same time windows)', ...
       'Flexible: bond turnover high while membership turnover low'})
legend('Location', 'best')
grid on

%% ---------------- VISUAL 2: PAIRED BARCODE ----------------
figure('Name', sprintf('Paired barcode at lambda_cross = %.3g', barcode_lambda), ...
       'Position', [100 100 1400 700])

% Top: membership barcode (which beads are in the persistent >=3 cluster)
subplot(2,1,1)
imagesc(barcode_trace.times, 1:barcode_trace.N, double(barcode_trace.membership_mat))
colormap(gca, [1 1 1; 0 0.4 0.8])
set(gca, 'YDir', 'reverse')
ylabel('bead #')
title(sprintf(['Membership: which beads are in the persistent (\\geq%d) cluster' ...
               '   (\\lambda_{cross} = %.3g)'], min_cluster_size, barcode_lambda))

% Bottom: bond barcode (which pairs among clustered beads are bonded)
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
title('Bonds: which specific pairs are crosslinked (flexible = stable band above, flickering below)')

%% ---------------- VISUAL 3: TURNOVER QUADRANT SCATTER ----------------
figure('Name', 'Regime quadrant', 'Position', [150 150 800 650])
scatter(bond_turnover, membership_turnover, 120, log10(lambda_cross_list), 'filled')
hold on
for s = 1:n_sweep
    text(bond_turnover(s)*1.05, membership_turnover(s), ...
         sprintf('%.2g', lambda_cross_list(s)), 'FontSize', 9);
end
xlabel('Bond turnover (fraction of bonds changing / s)')
ylabel('Membership turnover (fraction of cluster membership changing / s)')
title({'Regime quadrant', ...
       'RIGID: low-low   |   FLEXIBLE: high bond turnover, low membership turnover   |   AMORPHIC: high-high'})
cb = colorbar; cb.Label.String = 'log_{10}(\lambda_{cross})';
grid on

%% ---------------- SUMMARY TABLE ----------------
fprintf('\n========== REGIME SUMMARY ==========\n');
fprintf('%-12s %-16s %-14s %-18s %-16s %-14s\n', ...
    'lambda', 'memb_persist(s)', 'bond_persist(s)', 'memb_turnover(/s)', 'bond_turnover(/s)', 'frac_clustered');
fprintf('%s\n', repmat('-', 1, 95));
for s = 1:n_sweep
    fprintf('%-12.3g %-16.3f %-14.3f %-18.3f %-16.3f %-14.3f\n', ...
        lambda_cross_list(s), membership_persist(s), bond_persist(s), ...
        membership_turnover(s), bond_turnover(s), frac_time_clustered(s));
end
fprintf('\nInterpretation guide:\n');
fprintf('  RIGID    : memb_persist high, bond_persist high (ratio ~ 1)\n');
fprintf('  FLEXIBLE : memb_persist high, bond_persist low  (ratio >> 1)\n');
fprintf('  AMORPHIC : memb_persist low  (no stable >=%d cluster)\n', min_cluster_size);
fprintf('\nThe ratio memb_persist / bond_persist is the key flexible-regime signal:\n');
ratio = membership_persist ./ max(bond_persist, eps);
for s = 1:n_sweep
    fprintf('  lambda = %-8.3g  memb/bond ratio = %.2f\n', lambda_cross_list(s), ratio(s));
end

save('membership_bond_regime.mat', 'lambda_cross_list', 'membership_persist', ...
     'bond_persist', 'membership_turnover', 'bond_turnover', 'frac_time_clustered', ...
     'barcode_trace', 'barcode_lambda', 'sigma', 'k_cross', 'r_thresh', ...
     'min_cluster_size', 'fine_stride_s');
fprintf('\nSaved to membership_bond_regime.mat\n');


%% =========================================================================
function trace = build_membership_bond_trace(r, r_thresh, min_size, fine_stride_s)
% Walk the trajectory at a fine stride. At each sample:
%   - find the largest proximity cluster of size >= min_size (membership)
%   - record which of the bondable pairs are currently bonded (bond network)
% Then compute:
%   - membership persistence: how long the SAME bead-set (Jaccard>0.5) lasts
%   - bond persistence:       mean lifetime of an individual bond, from B
%   - membership/bond turnover rates (fraction changing per second)
%   - fraction of time a >= min_size cluster exists at all

trace = struct();
trace.N = r.N;

N      = r.N;
nsteps = size(r.traj, 3);
stride = max(1, round(fine_stride_s / r.dt));
idx    = stride : stride : nsteps;
nfine  = length(idx);
dt_fine = stride * r.dt;

membership_mat = false(N, nfine);     % bead in persistent cluster?
bond_mat       = false(N, N, nfine);  % pair bonded?
sets           = cell(1, nfine);      % membership bead-set at each sample

% Bonds: prefer fine-grained storage (B_fine) if available, aligned to the
% same physical times as the membership samples. Fall back to coarse
% B_history (5 s) only if fine storage is absent.
use_fine_bonds = isfield(r, 'B_fine') && ~isempty(r.B_fine);
if ~use_fine_bonds
    coarse_store = round(r.snap_times(1) / r.dt);
end

for k = 1:nfine
    t_idx  = idx(k);
    x_snap = r.traj(:, :, t_idx);
    d_mat  = sqrt((x_snap(:,1) - x_snap(:,1)').^2 + (x_snap(:,2) - x_snap(:,2)').^2);
    A      = (d_mat < r_thresh) & (d_mat > 0);
    G      = graph(A);
    bins   = conncomp(G);
    sizes  = accumarray(bins(:), 1);

    % largest cluster of size >= min_size (if any)
    big = find(sizes >= min_size);
    if isempty(big)
        sets{k} = [];
    else
        [~, which] = max(sizes(big));
        lbl = big(which);
        members = find(bins(:) == lbl);
        sets{k} = members;
        membership_mat(members, k) = true;
    end

    % bond network at this time
    if use_fine_bonds
        % B_fine is stored every (B_fine_times(1)/dt) steps; map t_idx to it
        bstride = round(r.B_fine_times(1) / r.dt);
        bidx    = min(max(1, round(t_idx / bstride)), size(r.B_fine, 3));
        bond_mat(:, :, k) = r.B_fine(:, :, bidx);
    else
        cidx = min(max(1, round(t_idx / coarse_store)), length(r.B_history));
        if ~isempty(r.B_history{cidx})
            bond_mat(:, :, k) = (r.B_history{cidx} > 0);
        end
    end
end

trace.times          = idx * r.dt;
trace.membership_mat = membership_mat;
trace.bond_mat       = bond_mat;

% --- membership persistence: run-length of stable >=3 bead-set ---
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
    trace.membership_persist = 0;
else
    trace.membership_persist = median(run_lengths) * dt_fine;
end

% --- bond persistence: mean individual bond lifetime (from recorded data) ---
if ~isempty(r.bond_lifetimes)
    trace.bond_persist = mean(r.bond_lifetimes);
else
    trace.bond_persist = 0;
end

% --- turnover rates (fraction changing per second) ---
% membership turnover: avg Jaccard distance between consecutive samples / dt
memb_changes = [];
for k = 2:nfine
    if isempty(sets{k}) && isempty(sets{k-1}), continue; end
    inter = intersect(sets{k}, sets{k-1});
    uni   = union(sets{k}, sets{k-1});
    jdist = 1 - numel(inter)/max(numel(uni),1);
    memb_changes(end+1) = jdist; %#ok<AGROW>
end
if isempty(memb_changes)
    trace.membership_turnover = 0;
else
    trace.membership_turnover = mean(memb_changes) / dt_fine;
end

% bond turnover: avg fraction of the bond set that changes between samples / dt
bond_changes = [];
for k = 2:nfine
    b1 = squeeze(bond_mat(:,:,k-1));
    b2 = squeeze(bond_mat(:,:,k));
    n1 = nnz(triu(b1,1));
    n2 = nnz(triu(b2,1));
    if n1 == 0 && n2 == 0, continue; end
    changed = nnz(triu(xor(b1,b2),1));
    bond_changes(end+1) = changed / max(n1+n2, 1); %#ok<AGROW>
end
if isempty(bond_changes)
    trace.bond_turnover = 0;
else
    trace.bond_turnover = mean(bond_changes) / dt_fine;
end

% --- fraction of time a >= min_size cluster exists ---
trace.frac_time_clustered = mean(cellfun(@(s) ~isempty(s), sets));
end