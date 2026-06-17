%% SINGLE-RUN DRIVER: one tau_cross value with full diagnostic plots
%  Use this for exploratory runs / debugging. For timescale analysis,
%  use sweep_tau_cross.m instead.

clear; clc; close all;

tau_cross    = 1;
steps        = 200000;
warmup_steps = 20000;

opts = struct('save_traj', true);
out  = run_hult_2d(tau_cross, steps, warmup_steps, opts);

% Unpack
N               = out.N;
R_nuc           = 175;
traj            = out.traj;
B_history       = out.B_history;
bond_fraction   = out.bond_fraction;
bond_lifetimes  = out.bond_lifetimes;
largest_cluster = out.largest_cluster;
cluster_sizes   = out.cluster_sizes;
snap_times      = out.snap_times;
x_final         = out.x_final;

%% ---------------- PLOTS ----------------

% --- Initial vs final chain configuration ---
figure
subplot(1,2,1); hold on
theta_disk = linspace(0, 2*pi, 200);
plot(R_nuc*cos(theta_disk), R_nuc*sin(theta_disk), 'k-')
plot(traj(:,1,1), traj(:,2,1), 'b-o', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
plot(traj(1,1,1), traj(1,2,1), 'gs', 'MarkerSize', 12, 'LineWidth', 2)
plot(traj(N,1,1), traj(N,2,1), 'rs', 'MarkerSize', 12, 'LineWidth', 2)
axis equal; grid on
xlabel('x (nm)'); ylabel('y (nm)')
title('Post-warmup initial configuration')

subplot(1,2,2); hold on
plot(R_nuc*cos(theta_disk), R_nuc*sin(theta_disk), 'k-')
plot(x_final(:,1), x_final(:,2), 'b-o', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
B_final = out.B_final;
for i = 1:N
    for j = i+1:N
        if B_final(i,j) == 1
            plot([x_final(i,1), x_final(j,1)], [x_final(i,2), x_final(j,2)], ...
                 'r-', 'LineWidth', 1.5)
        end
    end
end
plot(x_final(1,1), x_final(1,2), 'gs', 'MarkerSize', 12, 'LineWidth', 2)
plot(x_final(N,1), x_final(N,2), 'rs', 'MarkerSize', 12, 'LineWidth', 2)
axis equal; grid on
xlabel('x (nm)'); ylabel('y (nm)')
title('Final configuration (red = active crosslinks)')

% --- Time series of Rg and n_bonds ---
figure
subplot(2,1,1)
plot(out.time, out.Rg_t, 'b'); grid on
xlabel('Time (s)'); ylabel('R_g (nm)')
title(sprintf('R_g(t)  |  \\tau_{cross} = %.3g', tau_cross))
subplot(2,1,2)
plot(out.time, out.n_bonds_t, 'r'); grid on
xlabel('Time (s)'); ylabel('Active bonds')
title('Number of active crosslinks')

% --- Cluster growth ---
figure
plot(snap_times, largest_cluster, 'b', 'LineWidth', 2)
xlabel('Time (s)'); ylabel('Largest cluster size')
title(sprintf('Cluster growth (Hult hard-cutoff, 2D, \\tau_{cross}=%.3g)', tau_cross))
grid on

% --- Cluster size distribution ---
figure
histogram(cluster_sizes(cluster_sizes > 0), 'BinEdges', 0.5:1:N+0.5, 'FaceColor', 'b')
xlabel('Cluster size'); ylabel('Frequency')
title('Cluster size distribution'); grid on

% --- Bond lifetime distribution ---
figure
if ~isempty(bond_lifetimes)
    histogram(bond_lifetimes, 30, 'FaceColor', 'b')
    xline(mean(bond_lifetimes), '--k', 'LineWidth', 1.5)
    xlabel('Bond lifetime (s)'); ylabel('Count')
    title(sprintf('Bond lifetimes  |  mean = %.3f s, median = %.3f s', ...
        mean(bond_lifetimes), median(bond_lifetimes)))
else
    title('No bonds formed')
end
grid on

% --- Bond fraction heatmap ---
figure
bf = bond_fraction;
mask = abs((1:N)' - (1:N)) <= 1;
bf(mask) = NaN;
imagesc(bf); colorbar; colormap('hot'); axis square
xlabel('Bead j'); ylabel('Bead i')
title('Fraction of time bonded per bead pair')

%% ---------------- SUMMARY ----------------
fprintf('\n========== 2D HULT HARD-CUTOFF SUMMARY ==========\n')
if ~isempty(bond_lifetimes)
    eligible = [];
    for i = 1:N
        for j = i+1:N
            if abs(i-j) > 1
                eligible(end+1) = bond_fraction(i,j); %#ok<AGROW>
            end
        end
    end
    fprintf('  N beads, R_nuc     : %d, %.0f nm\n', N, R_nuc)
    fprintf('  tau_cross          : %.3g\n', tau_cross)
    fprintf('  warmup steps       : %d\n', warmup_steps)
    fprintf('  production steps   : %d\n', steps)
    fprintf('  kon_eff, koff_eff  : %.3g, %.3g\n', out.kon_eff, out.koff_eff)
    fprintf('  Total bonds formed : %d\n',     length(bond_lifetimes))
    fprintf('  Mean lifetime      : %.4f s\n', mean(bond_lifetimes))
    fprintf('  Median lifetime    : %.4f s\n', median(bond_lifetimes))
    fprintf('  Max lifetime       : %.4f s\n', max(bond_lifetimes))
    fprintf('  Std deviation      : %.4f s\n', std(bond_lifetimes))
    fprintf('  Mean bond fraction : %.4f\n',   mean(eligible))
    fprintf('  Max bond fraction  : %.4f\n',   max(eligible))
    fprintf('  Mean R_g           : %.2f nm\n', mean(out.Rg_t))
    fprintf('  Mean N_bonds       : %.2f\n',   mean(out.n_bonds_t))
else
    fprintf('  No bonds recorded\n')
end
fprintf('=================================================\n')