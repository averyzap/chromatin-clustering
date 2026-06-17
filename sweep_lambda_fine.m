%% DRIVER: fine-grained lambda_cross sweep + regime-distinguishing diagnostics
%
%  Two goals:
%    1. Increase resolution across the regime crossover by adding intermediate
%       lambda_cross points (0.3 and 3) to the existing 5-point sweep.
%    2. Replot the key observables against the dimensionless ratio
%       tau(N_bonds) / tau(R_g) instead of lambda_cross, so the regime
%       boundary is anchored to the bond-vs-polymer timescale ratio
%       (the physically meaningful quantity, connecting directly to Anna's
%       alpha parameter).
%
%  Adds a THIRD independent regime signature beyond cluster size and
%  autocorrelation crossover: cluster-identity persistence -- how long the
%  same set of beads stays mutually clustered. This should be long in the
%  rigid regime, short in the flexible regime.

clear; clc; close all;

%% ---------------- SWEEP CONFIGURATION ----------------
lambda_cross_list = [0.01, 0.1, 0.3, 1, 3, 10, 100];
steps_base        = 200000;
warmup_steps      = 20000;
sigma             = 60;
k_cross           = 0.1;

steps_list = round(steps_base * max(1, sqrt(lambda_cross_list)));

n_sweep = length(lambda_cross_list);
results = cell(1, n_sweep);

%% ---------------- (1) REFERENCE RUN: NO CROSSLINKS ----------------
fprintf('\n========== REFERENCE RUN (k_on = 0) ==========\n');
ref_opts = struct('k_on', 0, 'k_off', 0.5, 'sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 42, 'save_traj', false);
ref = run_hult_2d(1, steps_base, warmup_steps, ref_opts);

[tau_p, acf_ref, lags_ref] = autocorr_time(ref.Rg_t, ref.dt);
fprintf('Polymer relaxation time tau_p = %.3f s\n', tau_p);

%% ---------------- (2+3) SWEEP + PER-RUN ANALYSIS ----------------
% We analyze each run immediately after it finishes, while its (large)
% trajectory is still in memory, then strip the trajectory before storing.
% This keeps peak memory to one trajectory at a time and keeps the saved
% .mat small. Persistence MUST be computed here because it needs r.traj.
fprintf('\n========== lambda_cross FINE SWEEP ==========\n');

mean_lifetime      = zeros(1, n_sweep);
mean_n_bonds       = zeros(1, n_sweep);
mean_Rg            = zeros(1, n_sweep);
mean_largest_cl    = zeros(1, n_sweep);
tau_Rg             = zeros(1, n_sweep);
tau_nbonds         = zeros(1, n_sweep);
cluster_id_persist = zeros(1, n_sweep);

for s = 1:n_sweep
    lc = lambda_cross_list(s);
    fprintf('\n--- Sweep %d/%d: lambda_cross = %.3g ---\n', s, n_sweep, lc);
    opts = struct('sigma', sigma, 'k_cross', k_cross, ...
                  'seed', 200 + s, 'save_traj', true);
    r = run_hult_2d(lc, steps_list(s), warmup_steps, opts);

    % analyze while trajectory is in hand
    if ~isempty(r.bond_lifetimes)
        mean_lifetime(s) = mean(r.bond_lifetimes);
    end
    mean_n_bonds(s)    = mean(r.n_bonds_t);
    mean_Rg(s)         = mean(r.Rg_t);
    mean_largest_cl(s) = mean(r.largest_cluster);
    tau_Rg(s)          = autocorr_time(r.Rg_t,      r.dt);
    tau_nbonds(s)      = autocorr_time(r.n_bonds_t, r.dt);
    cluster_id_persist(s) = compute_cluster_identity_persistence(r);

    % strip the bulky trajectory before storing (persistence already computed)
    r.traj = [];
    results{s} = r;
end

tau_ratio = tau_nbonds ./ tau_Rg;

%% ---------------- (4) PLOTS ----------------

% --- Standard: observables vs lambda_cross ---
figure('Name', 'Fine sweep: observables vs lambda_cross', 'Position', [50 50 1200 800])

subplot(2,3,1)
loglog(lambda_cross_list, mean_lifetime, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on
loglog(lambda_cross_list, lambda_cross_list ./ 0.5, 'k--', 'LineWidth', 1)
xlabel('\lambda_{cross}'); ylabel('Mean bond lifetime (s)')
legend('measured', '1/k_{off,eff}', 'Location', 'best')
title('Bond lifetime'); grid on

subplot(2,3,2)
semilogx(lambda_cross_list, mean_n_bonds, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('\lambda_{cross}'); ylabel('\langle N_{bonds} \rangle')
title('Mean active bonds'); grid on

subplot(2,3,3)
semilogx(lambda_cross_list, mean_largest_cl, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('\lambda_{cross}'); ylabel('\langle largest cluster \rangle')
title('Largest cluster size'); grid on

subplot(2,3,4)
loglog(lambda_cross_list, tau_Rg, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on
loglog(lambda_cross_list, tau_nbonds, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r')
yline(tau_p, '--k', sprintf('\\tau_p = %.2f s', tau_p))
xlabel('\lambda_{cross}'); ylabel('Autocorrelation time (s)')
legend('\tau(R_g)', '\tau(N_{bonds})', 'Location', 'best')
title('Observable autocorrelation times'); grid on

subplot(2,3,5)
semilogx(lambda_cross_list, mean_Rg, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on
yline(mean(ref.Rg_t), '--k', 'no crosslinks')
xlabel('\lambda_{cross}'); ylabel('\langle R_g \rangle (nm)')
title('Chain compaction'); grid on

subplot(2,3,6)
semilogx(lambda_cross_list, cluster_id_persist, 'mo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'm')
xlabel('\lambda_{cross}'); ylabel('Cluster identity persistence (s)')
title('How long same bead-set stays clustered'); grid on

% --- Observables vs dimensionless ratio ---
figure('Name', 'Regime diagnostics: observables vs tau(N_bonds)/tau(R_g)', ...
       'Position', [100 100 1400 500])

subplot(1,3,1)
semilogx(tau_ratio, mean_largest_cl, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b')
hold on
xline(1, '--k', '\tau_{bonds} = \tau_p')
xlabel('\tau(N_{bonds}) / \tau(R_g)')
ylabel('\langle largest cluster \rangle')
title('Cluster size vs dimensionless timescale ratio')
grid on

subplot(1,3,2)
semilogx(tau_ratio, mean_n_bonds, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b')
hold on
xline(1, '--k', '\tau_{bonds} = \tau_p')
xlabel('\tau(N_{bonds}) / \tau(R_g)')
ylabel('\langle N_{bonds} \rangle')
title('Mean bonds vs dimensionless timescale ratio')
grid on

subplot(1,3,3)
semilogx(tau_ratio, cluster_id_persist, 'mo-', 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold on
xline(1, '--k', '\tau_{bonds} = \tau_p')
xlabel('\tau(N_{bonds}) / \tau(R_g)')
ylabel('Cluster identity persistence (s)')
title('Identity persistence vs dimensionless timescale ratio')
grid on

% --- Headline three-signature plot ---
figure('Name', 'Three independent regime signatures', 'Position', [150 150 900 600])
norm_cluster  = (mean_largest_cl - min(mean_largest_cl)) / (max(mean_largest_cl) - min(mean_largest_cl) + eps);
norm_persist  = (cluster_id_persist - min(cluster_id_persist)) / (max(cluster_id_persist) - min(cluster_id_persist) + eps);
norm_taubonds = log10(tau_nbonds / tau_p);
norm_taubonds = (norm_taubonds - min(norm_taubonds)) / (max(norm_taubonds) - min(norm_taubonds) + eps);

semilogx(lambda_cross_list, norm_cluster,  'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Largest cluster size'); hold on
semilogx(lambda_cross_list, norm_persist,  'mo-', 'LineWidth', 2, 'MarkerFaceColor', 'm', 'DisplayName', 'Cluster identity persistence');
semilogx(lambda_cross_list, norm_taubonds, 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'log[\tau(N_{bonds})/\tau_p]');
xlabel('\lambda_{cross}')
ylabel('Normalized signal (0 = min, 1 = max)')
title('Three independent regime signatures on one axis')
legend('Location', 'best')
grid on

%% ---------------- SUMMARY ----------------
fprintf('\n========== FINE-SWEEP SUMMARY ==========\n');
fprintf('Polymer relaxation time tau_p = %.3f s\n', tau_p);
fprintf('sigma = %.0f nm, k_cross = %.3g\n\n', sigma, k_cross);
fprintf('%-12s %-12s %-14s %-12s %-12s %-12s %-12s\n', ...
    'lambda_cross', '<largest>', '<bond life>', 'tau_ratio', '<N_bonds>', '<R_g>', 'persist(s)');
fprintf('%s\n', repmat('-', 1, 100));
for s = 1:n_sweep
    fprintf('%-12.3g %-12.2f %-14.3f %-12.3f %-12.2f %-12.1f %-12.3f\n', ...
        lambda_cross_list(s), mean_largest_cl(s), mean_lifetime(s), ...
        tau_ratio(s), mean_n_bonds(s), mean_Rg(s), cluster_id_persist(s));
end

save('lambda_fine_sweep_results.mat', 'results', 'lambda_cross_list', 'tau_p', ...
     'mean_lifetime', 'mean_n_bonds', 'mean_Rg', 'mean_largest_cl', ...
     'tau_Rg', 'tau_nbonds', 'tau_ratio', 'cluster_id_persist', ...
     'sigma', 'k_cross', 'ref');
fprintf('\nSaved results to lambda_fine_sweep_results.mat\n');


%% =========================================================================
function [tau, acf, lags] = autocorr_time(y, dt)
y = y(:) - mean(y);
n = length(y);
nfft = 2^nextpow2(2*n - 1);
Y    = fft(y, nfft);
acf_full = real(ifft(Y .* conj(Y)));
acf  = acf_full(1:n) ./ (n - (0:n-1).');
acf  = acf / acf(1);
lags = (0:n-1).' * dt;

idx = find(acf < 1/exp(1), 1, 'first');
if isempty(idx) || idx == 1
    tau = NaN; return
end
y1 = acf(idx-1); y2 = acf(idx);
frac = (1/exp(1) - y1) / (y2 - y1);
tau = lags(idx-1) + frac * dt;
end


%% =========================================================================
function persist_time = compute_cluster_identity_persistence(r)
% Cluster-identity persistence: at each snapshot, identify the largest
% spatial cluster (proximity-based, threshold = 25 nm) and record its
% member bead-set. Then count how many consecutive snapshots that set
% remains "substantially the same" (Jaccard >= 0.5). The persistence
% time is the median run-length in seconds.

persist_time = NaN;
if isempty(r.traj), return; end

N       = r.N;
nsteps  = size(r.traj, 3);

% Walk the trajectory at a FINE stride, independent of the coarse B_history
% storage interval. The coarse interval (5000 steps = 5 s) is far too long
% to resolve cluster-identity turnover, which happens on ~tau_p ~ 1 s. We
% recompute the largest-cluster membership directly from bead positions at
% every fine_stride steps. Default targets ~0.01 s resolution so the fast
% (flexible/frozen) end is not quantized by the sampling grid.
fine_stride = max(1, round(0.01 / r.dt));   % ~0.01 s between samples
snap_idx    = fine_stride : fine_stride : nsteps;
nsnaps      = length(snap_idx);
if nsnaps < 3, return; end

r_thresh = 25;

sets = cell(1, nsnaps);
for k = 1:nsnaps
    x_snap = r.traj(:, :, snap_idx(k));
    d_mat  = sqrt((x_snap(:,1) - x_snap(:,1)').^2 + (x_snap(:,2) - x_snap(:,2)').^2);
    A      = (d_mat < r_thresh) & (d_mat > 0);
    G      = graph(A);
    bins   = conncomp(G);
    sizes  = accumarray(bins(:), 1);
    [~, biggest] = max(sizes);
    sets{k} = find(bins == biggest);
end

run_lengths = [];
run_start = 1;
for k = 2:nsnaps
    if isempty(sets{k}) || isempty(sets{k-1})
        run_lengths(end+1) = (k - 1) - run_start + 1; %#ok<AGROW>
        run_start = k;
        continue;
    end
    inter = intersect(sets{k}, sets{k-1});
    uni   = union(sets{k}, sets{k-1});
    jac   = length(inter) / max(length(uni), 1);
    if jac < 0.5
        run_lengths(end+1) = (k - 1) - run_start + 1; %#ok<AGROW>
        run_start = k;
    end
end
run_lengths(end+1) = nsnaps - run_start + 1;

snap_dt = fine_stride * r.dt;
persist_time = median(run_lengths) * snap_dt;
end