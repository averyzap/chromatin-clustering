function out = run_hult_2d(lambda_cross, steps, warmup_steps, opts)
% RUN_HULT_2D  Tethered polymer simulation (2D) with Gaussian binding kernel.
%   (Filename kept for repo continuity; kernel is k_on(r) = k_on*exp(-r^2/sigma^2).)
%
%   out = run_hult_2d(lambda_cross, steps, warmup_steps)
%   out = run_hult_2d(lambda_cross, steps, warmup_steps, opts)
%
%   Runs a warmup phase (no data collection) followed by a production phase.
%
%   lambda_cross is the dimensionless rate-rescaling factor.
%     k_on_eff  = k_on  / lambda_cross
%     k_off_eff = k_off / lambda_cross
%     lambda_cross < 1 -> fast bonds (rigid limit)
%     lambda_cross = 1 -> baseline
%     lambda_cross > 1 -> slow bonds (amorphic limit)
%
%   dt is auto-tightened so max(k_on_eff, k_off_eff)*dt stays below 0.05.
%   When dt is reduced, steps and warmup_steps are scaled up so the total
%   simulated time is preserved.
%
%   opts (optional struct) fields:
%     .N             bead count (default 6)
%     .R_nuc         disk radius, nm (default 175)
%     .sigma         Gaussian kernel width, nm (default 60)
%     .dt            requested timestep (default 1e-3; may be auto-reduced)
%     .k_on, .k_off  base rates (defaults 1.0, 0.5)
%     .k_cross       crosslink spring stiffness (default 0.1)
%     .seed          RNG seed (default: shuffle)
%     .verbose       print progress (default true)
%     .save_traj     keep full trajectory (default true)
%     .bond_stride   if > 0, also store the bond matrix B every bond_stride
%                    steps into out.B_fine (with out.B_fine_times). Use this
%                    for fine-grained bond barcodes. Default 0 (off).

if nargin < 4, opts = struct(); end
if ~isfield(opts, 'N'),           opts.N           = 6;     end
if ~isfield(opts, 'R_nuc'),       opts.R_nuc       = 175;   end
if ~isfield(opts, 'sigma'),       opts.sigma       = 60;    end
if ~isfield(opts, 'dt'),          opts.dt          = 1e-3;  end
if ~isfield(opts, 'k_on'),        opts.k_on        = 1.0;   end
if ~isfield(opts, 'k_off'),       opts.k_off       = 0.5;   end
if ~isfield(opts, 'k_cross'),     opts.k_cross     = 0.1;   end
if ~isfield(opts, 'verbose'),     opts.verbose     = true;  end
if ~isfield(opts, 'save_traj'),   opts.save_traj   = true;  end
if ~isfield(opts, 'bond_stride'), opts.bond_stride = 0;     end
if isfield(opts, 'seed'),         rng(opts.seed); else, rng('shuffle'); end

%% ---------------- PARAMETERS ----------------
N     = opts.N;
R_nuc = opts.R_nuc;
dt_requested = opts.dt;

zeta = 2.5e-3;
kBT  = 4.1;

% chain WLC
Lp    = 50;
Nk    = 17;
R0    = Nk*(2*Lp);
alpha = 0.2176;

% Gaussian binding kernel: k_on(r) = k_on * exp(-r^2 / sigma^2)
sigma = opts.sigma;

% Numerical pruning cutoff: pairs farther than r_eligible are skipped because
% the Gaussian is negligible there. Bonded pairs may drift past it without
% breaking; it only gates bond formation.
r_eligible = 90;

cEV = 8.305e-5;
aEV = 3.268e-5;

k_on    = opts.k_on;
k_off   = opts.k_off;
k_cross = opts.k_cross;

k_on_eff  = k_on  / lambda_cross;
k_off_eff = k_off / lambda_cross;

% Auto-scale dt to keep the Poisson rate approximation valid.
p_target = 0.05;
rate_max = max(k_on_eff, k_off_eff);
if rate_max > 0
    dt_max = p_target / rate_max;
else
    dt_max = Inf;
end

if dt_requested > dt_max
    dt = dt_max;
    scale        = dt_requested / dt;
    steps        = round(steps * scale);
    warmup_steps = round(warmup_steps * scale);
    if opts.verbose
        fprintf(['[auto-dt] lambda_cross=%.3g forces dt=%.2e (requested ' ...
                 '%.2e); steps rescaled to %d (warmup %d) to preserve T_total.\n'], ...
                lambda_cross, dt, dt_requested, steps, warmup_steps);
    end
else
    dt = dt_requested;
end

noise = sqrt(2*kBT*dt/zeta);

% bond_stride is specified in STEPS at the requested dt; rescale if dt changed
bond_stride = opts.bond_stride;
if bond_stride > 0 && dt_requested > dt
    bond_stride = max(1, round(bond_stride * (dt_requested / dt)));
end

%% ---------------- INITIAL CONFIGURATION (SAW between pinned ends) ----------------
x = zeros(N, 2);
x(1, :) = [-R_nuc, 0];
x(N, :) = [ R_nuc, 0];

target_step = 70;
min_separation = 15;

max_attempts = 5000;
success = false;
for attempt = 1:max_attempts
    x_try = x;
    ok = true;
    for i = 2:N-1
        remaining = N - i + 1;
        target = (x(N,:) - x_try(i-1,:)) / remaining;
        target_dir = target / max(norm(target), eps);
        angle = atan2(target_dir(2), target_dir(1)) + 0.8 * (rand - 0.5) * pi;
        step = target_step * [cos(angle), sin(angle)];
        candidate = x_try(i-1,:) + step;
        if norm(candidate) > R_nuc - 10 || ...
                any(vecnorm(x_try(1:i-1,:) - candidate, 2, 2) < min_separation)
            ok = false;
            break
        end
        x_try(i,:) = candidate;
    end
    if ok
        if norm(x_try(N-1,:) - x(N,:)) < 3 * target_step
            x = x_try;
            success = true;
            break
        end
    end
end
if ~success
    warning('SAW init failed after %d attempts; using straight-line fallback', max_attempts);
    for i = 2:N-1
        x(i,:) = x(1,:) + (i-1)/(N-1) * (x(N,:) - x(1,:));
    end
end

%% ---------------- BOND BOOKKEEPING ----------------
B          = zeros(N);
bound      = zeros(N, 1);
bond_start = zeros(N);

%% ---------------- WARMUP LOOP (no data collection) ----------------
if opts.verbose
    fprintf('[lambda_cross=%.3g] Warmup: %d steps...\n', lambda_cross, warmup_steps);
end
tic
for t = 1:warmup_steps
    [x, B, bound, bond_start, ~] = step_dynamics( ...
        x, B, bound, bond_start, t, N, R_nuc, R0, alpha, k_cross, ...
        cEV, aEV, r_eligible, sigma, k_on_eff, k_off_eff, dt, zeta, noise, false);
end
if opts.verbose
    fprintf('[lambda_cross=%.3g] Warmup done (%.1f s wall time)\n', lambda_cross, toc);
end

% Reset bond_start so currently-active bonds start counting from production t=1
[ii, jj] = find(triu(B, 1));
bond_start = zeros(N);
for k = 1:length(ii)
    bond_start(ii(k), jj(k)) = 1;
end

%% ---------------- PRODUCTION LOOP ----------------
bond_lifetimes = [];
time_bonded    = zeros(N);
n_bonds_t      = zeros(steps, 1);
Rg_t           = zeros(steps, 1);
Ree_t          = zeros(steps, 1);

store_interval = 5000;
num_snaps      = floor(steps / store_interval);
B_history       = cell(1, num_snaps);
if opts.save_traj
    traj = zeros(N, 2, steps);
else
    traj = [];
end
largest_cluster = zeros(num_snaps, 1);
cluster_sizes   = zeros(num_snaps, N);
snap_idx = 1;

% optional fine-grained bond storage
if bond_stride > 0
    n_bond_fine = floor(steps / bond_stride);
    B_fine       = false(N, N, n_bond_fine);
    B_fine_times = (1:n_bond_fine) * bond_stride * dt;
    bond_fine_idx = 1;
else
    B_fine = [];
    B_fine_times = [];
end

if opts.verbose
    fprintf('[lambda_cross=%.3g] Production: %d steps...\n', lambda_cross, steps);
end
tic
for t = 1:steps
    [x, B, bound, bond_start, released] = step_dynamics( ...
        x, B, bound, bond_start, t, N, R_nuc, R0, alpha, k_cross, ...
        cEV, aEV, r_eligible, sigma, k_on_eff, k_off_eff, dt, zeta, noise, true);

    if ~isempty(released)
        bond_lifetimes = [bond_lifetimes, released]; %#ok<AGROW>
    end

    time_bonded = time_bonded + B;
    n_bonds_t(t) = sum(B(:)) / 2;
    com = mean(x, 1);
    Rg_t(t)  = sqrt(mean(sum((x - com).^2, 2)));
    Ree_t(t) = norm(x(N,:) - x(1,:));

    if opts.save_traj
        traj(:,:,t) = x;
    end

    if bond_stride > 0 && mod(t, bond_stride) == 0 && bond_fine_idx <= n_bond_fine
        B_fine(:, :, bond_fine_idx) = (B > 0);
        bond_fine_idx = bond_fine_idx + 1;
    end

    if mod(t, store_interval) == 0
        B_history{snap_idx} = B;
        if opts.save_traj
            x_snap = traj(:, :, t);
        else
            x_snap = x;
        end
        d_mat = sqrt((x_snap(:,1) - x_snap(:,1)').^2 + (x_snap(:,2) - x_snap(:,2)').^2);
        A     = (d_mat < 25) & (d_mat > 0);
        G     = graph(A);
        bins  = conncomp(G);
        sizes = accumarray(bins(:), 1).';
        cluster_sizes(snap_idx, 1:length(sizes)) = sizes;
        largest_cluster(snap_idx) = max(sizes);
        snap_idx = snap_idx + 1;
    end
end

% capture bonds still active at end of run
for i = 1:N
    for j = i+1:N
        if B(i,j) == 1
            lifetime = (steps - bond_start(i,j)) * dt;
            bond_lifetimes(end+1) = lifetime; %#ok<AGROW>
        end
    end
end
if opts.verbose
    fprintf('[lambda_cross=%.3g] Production done (%.1f s wall time)\n', lambda_cross, toc);
end

%% ---------------- PACKAGE OUTPUT ----------------
out = struct();
out.lambda_cross    = lambda_cross;
out.steps           = steps;
out.warmup_steps    = warmup_steps;
out.dt              = dt;
out.N               = N;
out.k_on_eff        = k_on_eff;
out.k_off_eff       = k_off_eff;
out.bond_lifetimes  = bond_lifetimes;
out.time_bonded     = time_bonded;
out.bond_fraction   = (time_bonded / 2) / steps;
out.n_bonds_t       = n_bonds_t;
out.Rg_t            = Rg_t;
out.Ree_t           = Ree_t;
out.largest_cluster = largest_cluster;
out.cluster_sizes   = cluster_sizes;
out.B_history       = B_history;
out.B_final         = B;
out.x_final         = x;
out.traj            = traj;
out.B_fine          = B_fine;
out.B_fine_times    = B_fine_times;
out.snap_times      = (1:num_snaps) * store_interval * dt;
out.time            = (1:steps) * dt;

end

%% =========================================================================
function [x, B, bound, bond_start, released_lifetimes] = step_dynamics( ...
    x, B, bound, bond_start, t, N, R_nuc, R0, alpha, k_cross, ...
    cEV, aEV, r_eligible, sigma, k_on_eff, k_off_eff, dt, zeta, noise, record_lifetimes)
% One Brownian-dynamics step with Gaussian-kernel crosslink kinetics.
% k_on(r) = k_on_eff * exp(-r^2 / sigma^2). Pairs farther than r_eligible
% are skipped. One bond per bead.

released_lifetimes = [];

% --- unbinding pass ---
for i = 1:N
    for j = i+1:N
        if B(i,j) == 1 && rand < k_off_eff * dt
            if record_lifetimes
                lifetime = (t - bond_start(i,j)) * dt;
                released_lifetimes(end+1) = lifetime; %#ok<AGROW>
            end
            B(i,j) = 0; B(j,i) = 0;
            bound(i) = 0; bound(j) = 0;
            bond_start(i,j) = 0;
        end
    end
end

% --- binding pass: per-pair Gaussian on-rate, no closest-pair priority ---
for i = 1:N
    for j = i+1:N
        if abs(i-j) > 1 && B(i,j) == 0 && bound(i) == 0 && bound(j) == 0
            d = norm(x(i,:) - x(j,:));
            if d < r_eligible
                k_on_r = k_on_eff * exp(-(d^2) / (sigma^2));
                if rand < k_on_r * dt
                    B(i,j) = 1; B(j,i) = 1;
                    bound(i) = 1; bound(j) = 1;
                    bond_start(i,j) = t;
                end
            end
        end
    end
end

% --- FORCES ---
F = zeros(N, 2);

% chain WLC
for i = 1:N-1
    rij  = x(i+1,:) - x(i,:);
    rabs = norm(rij);
    if rabs > 0 && rabs < R0
        f_mag = alpha * (-1 + 1/(1 - rabs/R0)^2 + 4*rabs/R0);
        f_vec = f_mag * rij / rabs;
        F(i,:)   = F(i,:)   + f_vec;
        F(i+1,:) = F(i+1,:) - f_vec;
    end
end

% crosslink linear springs
for d = 1:2
    dx_d   = x(:,d).' - x(:,d);
    Fcross = k_cross * B .* dx_d;
    F(:,d) = F(:,d) + sum(Fcross, 2);
end

% excluded volume (Gaussian repulsion)
r2 = (x(:,1) - x(:,1).').^2 + (x(:,2) - x(:,2).').^2;
for d = 1:2
    dx_d   = x(:,d) - x(:,d).';
    FEV    = cEV * dx_d .* exp(-aEV * r2);
    F(:,d) = F(:,d) + sum(FEV, 2);
end

% --- INTEGRATE ---
xi = noise * randn(N, 2);
x  = x + (dt/zeta)*F + xi;

% wall enforcement
r_norms = vecnorm(x, 2, 2);
outside = r_norms > R_nuc;
if any(outside)
    x(outside, :) = R_nuc * x(outside, :) ./ r_norms(outside);
end

% re-pin endpoints
x(1, :) = [-R_nuc, 0];
x(N, :) = [ R_nuc, 0];

end