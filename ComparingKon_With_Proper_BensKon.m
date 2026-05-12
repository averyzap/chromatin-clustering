%% TETHERED POLYMER SIMULATION (1D)
% WITH:
% - Single binding site per bead
% - Logistic Ben-style kon kernel
% - Matched kon scaling at r = 20 nm
% - Proximity-based clustering

clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N = 6;                 
L = 120;               
steps = 200000;
dt = 1e-3;

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
kon0  = 1.0;
koff  = 0.5;
k_cross = 0.01;

use_ben_kon = false;

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- BEN KERNEL MATCHING ----------------
r_match = 20;

% Gaussian reference (your baseline)
kG = kon0 * exp(-(r_match^2)/(20^2));  % sigma=20 as before

% Logistic kernel at match point
r_scaled_match = r_match / 45;
a_match = 2 / (1 + exp(20*(r_scaled_match - 0.75)));

% scale factor so models match at r = 20 nm
scale_factor = kG / (kon0 * a_match);

fprintf('Ben model scale factor = %f\n', scale_factor);

%% ---------------- VISUALIZE KON FUNCTIONS ----------------
r_vals = linspace(0,80,500);

k_gauss = kon0 * exp(-(r_vals.^2)/(20^2));
k_ben = zeros(size(r_vals));

for i = 1:length(r_vals)

    r_scaled = r_vals(i)/45;
    a = 2 / (1 + exp(20*(r_scaled - 0.75)));

    k_ben(i) = scale_factor * kon0 * a;
end

figure
plot(r_vals, k_gauss, 'b', 'LineWidth',2); hold on
plot(r_vals, k_ben, 'r', 'LineWidth',2)
xlabel('Distance (nm)')
ylabel('k_{on}')
legend('Gaussian','Ben logistic (scaled)')
title('Matched k_{on} Functions')
grid on

%% ---------------- INITIAL CONDITION ----------------
x = linspace(0,L,N)';
x(1) = 0;
x(N) = L;

traj = zeros(N,steps);

B = zeros(N);
bound = zeros(N,1);

store_interval = 5000;
B_history = {};
snap_idx = 1;

%% ---------------- SIMULATION ----------------
tic
for t = 1:steps

    %% ----- CROSSLINK DYNAMICS -----
    for i = 1:N
        for j = i+1:N
            
            if abs(i-j) > 1

                dist = x(i) - x(j);

                % -------- BEN LOGISTIC KON --------
                if use_ben_kon
                    r_scaled = abs(dist)/45;
                    a_ij = 2 / (1 + exp(20*(r_scaled - 0.75)));
                    kon_ij = scale_factor * kon0 * a_ij;
                else
                    kon_ij = kon0 * exp(-(dist^2)/(20^2));
                end

                % clamp
                kon_ij = max(0, kon_ij);
                kon_ij = min(kon_ij, 1/dt);

                % binding (valency = 1)
                if B(i,j) == 0 && bound(i)==0 && bound(j)==0
                    if rand < kon_ij * dt
                        B(i,j) = 1;
                        B(j,i) = 1;
                        bound(i) = 1;
                        bound(j) = 1;
                    end

                % unbinding
                elseif B(i,j) == 1
                    if rand < koff * dt
                        B(i,j) = 0;
                        B(j,i) = 0;
                        bound(i) = 0;
                        bound(j) = 0;
                    end
                end

            end
        end
    end

    %% ----- FORCES -----
    F = zeros(N,1);

    % WLC
    r = diff(x);
    rabs = abs(r);
    FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);

    F(2:N)   = F(2:N)   - sign(r).*FWLC;
    F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

    % Crosslinks
    dx = x' - x;
    Fcross = k_cross * B .* dx;
    F = F + sum(Fcross,2);

    % Excluded volume
    dx = x - x';
    FEV = cEV * dx .* exp(-aEV * dx.^2);
    F = F + sum(FEV,2);

    % Noise
    xi = noise * randn(N,1);

    % Update
    x = x + (dt/zeta)*F + xi;

    % tethering
    x(1) = 0;
    x(N) = L;

    traj(:,t) = x;

    if mod(t, store_interval) == 0
        B_history{snap_idx} = B;
        snap_idx = snap_idx + 1;
    end
end
toc

%% ---------------- PROXIMITY CLUSTER ANALYSIS ----------------
num_snaps = length(B_history);
cluster_sizes = zeros(num_snaps, N);
largest_cluster = zeros(num_snaps,1);

r_thresh = 20;

for k = 1:num_snaps
    
    x_snap = traj(:, k*store_interval);

    dx = abs(x_snap - x_snap');
    A = (dx < r_thresh) & (dx > 0);

    G = graph(A);
    bins = conncomp(G);

    sizes = zeros(1,max(bins));
    for c = 1:max(bins)
        sizes(c) = sum(bins == c);
    end

    cluster_sizes(k,1:length(sizes)) = sizes;
    largest_cluster(k) = max(sizes);
end

%% ---------------- VISUALIZATION ----------------

figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_history{k})
    axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Bond Network')

figure
plot(largest_cluster,'LineWidth',2)
xlabel('Snapshot')
ylabel('Largest cluster')
title('Cluster Growth')
grid on

figure
histogram(cluster_sizes(:).*(cluster_sizes(:)>0))
xlabel('Cluster size')
ylabel('Frequency')
title('Cluster Distribution')

%% ---------------- FINAL SNAPSHOT ----------------
figure
hold on
plot(x, zeros(size(x)), 'ko','MarkerFaceColor','k')

for i = 1:N
    for j = i+1:N
        if B(i,j) == 1
            plot([x(i), x(j)], [0,0], 'r-','LineWidth',2)
        end
    end
end

title('Final Configuration')
ylim([-1 1])