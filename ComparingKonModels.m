%% TETHERED POLYMER SIMULATION (1D)
% WITH:
% - Single binding site per bead
% - Switchable kon models
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
sigma = 20;
k_cross = 0.01;

use_ben_kon = true;   % toggle model

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- MATCHING KON AT r = 20 nm ----------------
r_match = 20;

% Gaussian value
kG = kon0 * exp(-(r_match^2)/(sigma^2));

% Ben raw value
r_tilde = min(r_match/45, 0.999);
kB_raw = 10.8794 * (-1 + 1/(1 - r_tilde)^2 + 4*r_tilde);

% Scale factor
scale_factor = kG / kB_raw;

fprintf('Ben model scale factor = %f\n', scale_factor);

%% ---------------- VISUALIZE KON FUNCTIONS ----------------
r_vals = linspace(0,80,500);

k_gauss = kon0 * exp(-(r_vals.^2)/(sigma^2));
k_ben = zeros(size(r_vals));

for i = 1:length(r_vals)
    r_tilde = min(r_vals(i)/45, 0.999);
    k_ben(i) = scale_factor * 10.8794 * (-1 + 1/(1 - r_tilde)^2 + 4*r_tilde);
end

k_ben = max(0, k_ben);

figure
plot(r_vals, k_gauss, 'b', 'LineWidth',2); hold on
plot(r_vals, k_ben, 'r', 'LineWidth',2)
xlabel('Distance (nm)')
ylabel('k_{on}')
legend('Gaussian','Scaled Ben model')
title('Matched k_{on} Functions (r = 20 nm)')
grid on

%% ---------------- INITIAL CONDITION ----------------
x = linspace(0,L,N)';
x(1) = 0;
x(N) = L;

traj = zeros(N,steps);

B = zeros(N);
bound = zeros(N,1);   % valency constraint

%% ---------------- STORAGE ----------------
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

                % --- SELECT KON MODEL ---
                if use_ben_kon
                    r_tilde = min(abs(dist)/45, 0.999);
                    kon_ij = scale_factor * 10.8794 * ...
                        (-1 + 1/(1 - r_tilde)^2 + 4*r_tilde);
                else
                    kon_ij = kon0 * exp(-(dist^2)/(sigma^2));
                end

                % --- CLAMP (CRITICAL) ---
                kon_ij = max(0, kon_ij);
                kon_ij = min(kon_ij, 1/dt);

                % --- BINDING (VALENCY = 1) ---
                if B(i,j) == 0 && bound(i)==0 && bound(j)==0
                    if rand < kon_ij * dt
                        B(i,j) = 1;
                        B(j,i) = 1;
                        bound(i) = 1;
                        bound(j) = 1;
                    end

                % --- UNBINDING ---
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

    % Crosslink
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

    % Tethers
    x(1) = 0;
    x(N) = L;

    traj(:,t) = x;

    % Store bond network
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

% Bonds
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_history{k})
    axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Bond Network (Valency = 1)')

% Proximity graphs
figure
for k = 1:num_snaps
    
    x_snap = traj(:, k*store_interval);
    dx = abs(x_snap - x_snap');
    A = (dx < r_thresh) & (dx > 0);

    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    G = graph(A);
    plot(G,'Layout','force')
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Proximity-Based Networks')

% Largest cluster
figure
plot(largest_cluster,'LineWidth',2)
xlabel('Snapshot index')
ylabel('Largest cluster size')
title('Cluster Growth (Proximity)')
grid on

% Distribution
all_clusters = cluster_sizes(:);
all_clusters = all_clusters(all_clusters > 0);

figure
histogram(all_clusters)
xlabel('Cluster size')
ylabel('Frequency')
title('Cluster Size Distribution')

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

title('Final Configuration with Bonds')
ylim([-1 1])