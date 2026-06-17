%% DRIVER: k_cross sensitivity sweep at fixed lambda_cross
%
%  Tests whether cluster physics is smoothly tunable with k_cross or whether
%  k_cross = 0.1 happens to sit near a sharp transition.
%
%  Observables tracked per k_cross:
%    - Largest cluster size              (the key cluster signal)
%    - <N_bonds>                         (mean number of active bonds)
%    - Mean bonded-pair separation       (the cleanest physical readout of
%                                         whether bonds are doing work or
%                                         are just nominally tracked)
%    - <R_g>                             (chain compaction)
%
%  Held fixed:
%    - lambda_cross  (set to 0.01, where clustering is strongest)
%    - sigma         (Gaussian kernel width, 60 nm)
%    - N, R_nuc      (chain geometry)
%
%  Requires run_hult_2d.m to accept opts.k_cross.

clear; clc; close all;

%% ---------------- SWEEP CONFIGURATION ----------------
k_cross_list = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0];
lambda_cross = 0.01;
sigma        = 60;

steps_base   = 200000;
warmup_steps = 20000;

n_sweep = length(k_cross_list);
results = cell(1, n_sweep);

%% ---------------- SWEEP ----------------
fprintf('\n========== k_cross SENSITIVITY SWEEP ==========\n');
fprintf('lambda_cross = %.3g (fixed), sigma = %.0f nm (fixed)\n\n', lambda_cross, sigma);

for s = 1:n_sweep
    kc = k_cross_list(s);
    fprintf('\n--- Sweep %d/%d: k_cross = %.3g ---\n', s, n_sweep, kc);

    opts = struct('sigma',     sigma, ...
                  'k_cross',   kc, ...
                  'seed',      100 + s, ...
                  'save_traj', true);
    results{s} = run_hult_2d(lambda_cross, steps_base, warmup_steps, opts);
end

%% ---------------- POST-HOC: MEAN BONDED-PAIR SEPARATION ----------------
mean_bond_sep    = zeros(1, n_sweep);
mean_largest_cl  = zeros(1, n_sweep);
mean_n_bonds     = zeros(1, n_sweep);
mean_Rg          = zeros(1, n_sweep);
mean_lifetime    = zeros(1, n_sweep);
n_obs_pairs      = zeros(1, n_sweep);

for s = 1:n_sweep
    r = results{s};

    seps_all = [];
    if ~isempty(r.B_history) && ~isempty(r.traj)
        store_interval = round(r.snap_times(1) / r.dt);
        for k = 1:length(r.B_history)
            Bk = r.B_history{k};
            if isempty(Bk), continue; end
            t_idx = k * store_interval;
            if t_idx > size(r.traj, 3), break; end
            x_snap = r.traj(:, :, t_idx);
            [ii, jj] = find(triu(Bk, 1) == 1);
            for m = 1:length(ii)
                d = norm(x_snap(ii(m),:) - x_snap(jj(m),:));
                seps_all(end+1) = d; %#ok<AGROW>
            end
        end
    end
    if ~isempty(seps_all)
        mean_bond_sep(s) = mean(seps_all);
        n_obs_pairs(s)   = length(seps_all);
    else
        mean_bond_sep(s) = NaN;
        n_obs_pairs(s)   = 0;
    end

    if ~isempty(r.largest_cluster)
        mean_largest_cl(s) = mean(r.largest_cluster);
    end
    mean_n_bonds(s) = mean(r.n_bonds_t);
    mean_Rg(s)      = mean(r.Rg_t);
    if ~isempty(r.bond_lifetimes)
        mean_lifetime(s) = mean(r.bond_lifetimes);
    end
end

%% ---------------- PLOTS ----------------
figure('Name', 'k_cross sensitivity', 'Position', [50 50 1400 800])

subplot(2,2,1)
semilogx(k_cross_list, mean_largest_cl, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('k_{cross}'); ylabel('\langle largest cluster \rangle')
title(sprintf('Largest cluster vs k_{cross}  (\\lambda_{cross}=%.3g)', lambda_cross))
grid on

subplot(2,2,2)
semilogx(k_cross_list, mean_n_bonds, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('k_{cross}'); ylabel('\langle N_{bonds} \rangle')
title('Mean active bonds vs k_{cross}'); grid on

subplot(2,2,3)
semilogx(k_cross_list, mean_bond_sep, 'ro-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r')
hold on
yline(sigma, '--k', sprintf('\\sigma = %.0f nm', sigma))
xlabel('k_{cross}'); ylabel('Mean bonded-pair separation (nm)')
title('Bond tightness vs k_{cross}'); grid on

subplot(2,2,4)
semilogx(k_cross_list, mean_Rg, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('k_{cross}'); ylabel('\langle R_g \rangle (nm)')
title('Chain compaction vs k_{cross}'); grid on

%% ---------------- SUMMARY TABLE ----------------
fprintf('\n========== k_cross SWEEP SUMMARY ==========\n');
fprintf('lambda_cross = %.3g, sigma = %.0f nm\n\n', lambda_cross, sigma);
fprintf('%-10s %-14s %-14s %-18s %-12s %-10s\n', ...
    'k_cross', '<largest cl>', '<N_bonds>', '<bonded sep> (nm)', '<R_g> (nm)', '# pairs');
fprintf('%s\n', repmat('-', 1, 85));
for s = 1:n_sweep
    fprintf('%-10.3g %-14.2f %-14.2f %-18.1f %-12.1f %-10d\n', ...
        k_cross_list(s), mean_largest_cl(s), mean_n_bonds(s), ...
        mean_bond_sep(s), mean_Rg(s), n_obs_pairs(s));
end

save('kcross_sweep_results.mat', 'results', 'k_cross_list', 'lambda_cross', ...
     'sigma', 'mean_largest_cl', 'mean_n_bonds', 'mean_bond_sep', ...
     'mean_Rg', 'mean_lifetime', 'n_obs_pairs');
fprintf('\nSaved results to kcross_sweep_results.mat\n');