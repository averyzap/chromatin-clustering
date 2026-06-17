%% DRIVER: lambda_cross sweep + timescale analysis for 2D Gaussian-kernel simulation

clear; clc; close all;

%% ---------------- SWEEP CONFIGURATION ----------------
lambda_cross_list = [0.01, 0.1, 1, 10, 100];
steps_base     = 200000;
warmup_steps   = 20000;

% For slow bonds need longer runs to get bond-lifetime statistics
steps_list = round(steps_base * max(1, sqrt(lambda_cross_list)));

n_sweep = length(lambda_cross_list);
results = cell(1, n_sweep);

%% ---------------- (1) REFERENCE RUN: NO CROSSLINKS ----------------
fprintf('\n========== REFERENCE RUN (k_on = 0) ==========\n');
ref_opts = struct('k_on', 0, 'k_off', 0.5, 'seed', 42, 'save_traj', false);
ref = run_hult_2d(1, steps_base, warmup_steps, ref_opts);

% Polymer relaxation time tau_p from Rg autocorrelation
[tau_p, acf_ref, lags_ref] = autocorr_time(ref.Rg_t, ref.dt);
fprintf('Polymer relaxation time tau_p = %.3f s\n', tau_p);

%% ---------------- (2) SWEEP lambda_cross ----------------
fprintf('\n========== lambda_cross SWEEP ==========\n');
for s = 1:n_sweep
    lc = lambda_cross_list(s);
    fprintf('\n--- Sweep %d/%d: lambda_cross = %.3g ---\n', s, n_sweep, lc);
    opts = struct('seed', 100 + s, 'save_traj', false);
    results{s} = run_hult_2d(lc, steps_list(s), warmup_steps, opts);
end

%% ---------------- (3) ANALYSIS PER lambda_cross ----------------
mean_lifetime    = zeros(1, n_sweep);
median_lifetime  = zeros(1, n_sweep);
mean_n_bonds     = zeros(1, n_sweep);
mean_Rg          = zeros(1, n_sweep);
mean_largest_cl  = zeros(1, n_sweep);
tau_Rg           = zeros(1, n_sweep);
tau_nbonds       = zeros(1, n_sweep);
exp_fit_R2       = zeros(1, n_sweep);   % goodness of exponential fit to lifetimes

for s = 1:n_sweep
    r = results{s};
    if ~isempty(r.bond_lifetimes)
        mean_lifetime(s)   = mean(r.bond_lifetimes);
        median_lifetime(s) = median(r.bond_lifetimes);
        exp_fit_R2(s)      = exp_survival_fit(r.bond_lifetimes);
    end
    mean_n_bonds(s)    = mean(r.n_bonds_t);
    mean_Rg(s)         = mean(r.Rg_t);
    mean_largest_cl(s) = mean(r.largest_cluster);
    tau_Rg(s)          = autocorr_time(r.Rg_t, r.dt);
    tau_nbonds(s)      = autocorr_time(r.n_bonds_t, r.dt);
end

% Expected bond lifetime from rate alone (no geometric constraints)
tau_off_naive = 1 ./ (0.5 ./ lambda_cross_list);   % = lambda_cross / 0.5 = 2 * lambda_cross

%% ---------------- (4) PLOTS ----------------

% --- Reference Rg autocorrelation ---
figure('Name', 'Reference polymer relaxation')
subplot(1,2,1)
plot(ref.time, ref.Rg_t, 'b'); grid on
xlabel('Time (s)'); ylabel('R_g (nm)')
title('R_g(t) with crosslinks off')
subplot(1,2,2)
plot(lags_ref, acf_ref, 'b', 'LineWidth', 1.5); hold on
yline(1/exp(1), '--k', '1/e')
xline(tau_p, '--r', sprintf('\\tau_p = %.2f s', tau_p))
xlim([0, min(10*tau_p, lags_ref(end))])
xlabel('Lag (s)'); ylabel('ACF of R_g')
title('Polymer relaxation autocorrelation')
grid on

% --- Observables vs lambda_cross ---
figure('Name', 'Observables vs lambda_cross', 'Position', [100 100 1200 700])

subplot(2,3,1)
loglog(lambda_cross_list, mean_lifetime, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on
loglog(lambda_cross_list, tau_off_naive, 'k--', 'LineWidth', 1)
xlabel('\lambda_{cross}'); ylabel('Mean bond lifetime (s)')
legend('measured', '1/k_{off,eff}', 'Location', 'best')
title('Bond lifetime'); grid on

subplot(2,3,2)
semilogx(lambda_cross_list, mean_n_bonds, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('\lambda_{cross}'); ylabel('\langle N_{bonds} \rangle')
title('Mean number of active bonds'); grid on

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
semilogx(lambda_cross_list, exp_fit_R2, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
xlabel('\lambda_{cross}'); ylabel('R^2 of exponential fit')
ylim([0 1.05])
title('Bond lifetime: Poisson-like?'); grid on

% --- Bond lifetime survival curves (semilog) ---
figure('Name', 'Bond lifetime survival functions')
cmap = lines(n_sweep);
hold on
for s = 1:n_sweep
    lt = results{s}.bond_lifetimes;
    if isempty(lt), continue; end
    lt_sorted = sort(lt);
    surv = 1 - (1:length(lt_sorted)) / length(lt_sorted);
    semilogy(lt_sorted, surv, '-', 'Color', cmap(s,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\lambda_{cross}=%.3g', lambda_cross_list(s)))
end
set(gca, 'YScale', 'log')
xlabel('Bond lifetime (s)')
ylabel('Survival probability P(T > t)')
title('Bond lifetime survival (straight line = exponential)')
legend('Location', 'best')
grid on

% --- Crossover regime indicator ---
figure('Name', 'Regime crossover')
loglog(lambda_cross_list, mean_lifetime / tau_p, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
hold on
yline(1, '--k', '\tau_{off} = \tau_p')
xlabel('\lambda_{cross}'); ylabel('\tau_{off} / \tau_p')
title('Timescale separation regime')
grid on
text(lambda_cross_list(1), 0.1, 'quasi-static bonds (fast)', 'FontSize', 10)
text(lambda_cross_list(end-1), 10, 'frozen bonds (slow)', 'FontSize', 10)

%% ---------------- SUMMARY TABLE ----------------
fprintf('\n========== SWEEP SUMMARY ==========\n');
fprintf('Polymer relaxation time tau_p = %.3f s\n\n', tau_p);
fprintf('%-12s %-14s %-14s %-12s %-12s %-14s %-12s\n', ...
    'lambda_cross', 'tau_off (meas)', 'tau_off/tau_p', '<N_bonds>', '<Rg>', 'tau(Rg)', 'exp R^2');
fprintf('%s\n', repmat('-', 1, 100));
for s = 1:n_sweep
    fprintf('%-12.3g %-14.3f %-14.3f %-12.2f %-12.1f %-14.3f %-12.3f\n', ...
        lambda_cross_list(s), mean_lifetime(s), mean_lifetime(s)/tau_p, ...
        mean_n_bonds(s), mean_Rg(s), tau_Rg(s), exp_fit_R2(s));
end

% Save sweep results
save('gaussian_sweep_results.mat', 'results', 'lambda_cross_list', 'tau_p', ...
     'mean_lifetime', 'mean_n_bonds', 'mean_Rg', 'mean_largest_cl', ...
     'tau_Rg', 'tau_nbonds', 'exp_fit_R2', 'ref');
fprintf('\nSaved sweep results to gaussian_sweep_results.mat\n');

%% =========================================================================
function [tau, acf, lags] = autocorr_time(y, dt)
% AUTOCORR_TIME  Integrated/exponential autocorrelation time of a time series.
%   Returns tau such that ACF(tau) = 1/e (linear interpolation between samples).
%   Also returns the normalized one-sided ACF and corresponding lag times.

y = y(:) - mean(y);
n = length(y);
% Use FFT-based autocorrelation for speed
nfft = 2^nextpow2(2*n - 1);
Y    = fft(y, nfft);
acf_full = real(ifft(Y .* conj(Y)));
acf  = acf_full(1:n) ./ (n - (0:n-1).');
acf  = acf / acf(1);
lags = (0:n-1).' * dt;

% Find first crossing of 1/e
idx = find(acf < 1/exp(1), 1, 'first');
if isempty(idx) || idx == 1
    tau = NaN;
    return
end
% Linear interpolation between idx-1 and idx
y1 = acf(idx-1); y2 = acf(idx);
frac = (1/exp(1) - y1) / (y2 - y1);
tau = lags(idx-1) + frac * dt;
end

%% =========================================================================
function R2 = exp_survival_fit(lifetimes)
% EXP_SURVIVAL_FIT  R^2 of linear fit to log(survival) vs lifetime.
%   R^2 ~ 1 means lifetime distribution is exponential (Poisson unbinding).
%   R^2 < 1 indicates deviations (e.g. geometric frustration, sub/super-exp tails).

if length(lifetimes) < 10
    R2 = NaN; return
end
lt = sort(lifetimes(:));
surv = 1 - (1:length(lt)).' / length(lt);
% drop the last point (survival = 0 -> log undefined)
lt = lt(1:end-1);
surv = surv(1:end-1);
% drop any zeros for safety
mask = surv > 0;
lt = lt(mask); surv = surv(mask);
if length(lt) < 5, R2 = NaN; return; end

log_s = log(surv);
p = polyfit(lt, log_s, 1);
log_s_fit = polyval(p, lt);
ss_res = sum((log_s - log_s_fit).^2);
ss_tot = sum((log_s - mean(log_s)).^2);
R2 = 1 - ss_res / ss_tot;
end