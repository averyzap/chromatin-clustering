%% TETHERED POLYMER SIMULATION (1D)
% WITH:
% - Single binding site per bead
% - Switchable kon models
% - Proximity-based network analysis

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

% --- SWITCH FOR KON MODEL ---
use_ben_kon = false;   % true = Ben's model, false = Gaussian (Caitlin's)

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONDITION ----------------
x = linspace(0,L,N)';
x(1) = 0;
x(N) = L;

traj = zeros(N,steps);

% Crosslink matrix
B = zeros(N);

% --- NEW: Occupancy tracker (valency = 1) ---
bound = zeros(N,1);

%% ---------------- NETWORK STORAGE ----------------
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

                % ----- KON MODEL -----
                if use_ben_kon
                    r_tilde = min(abs(dist)/45, 0.999);
                    kon_ij = 10.8794 * (-1 + 1/(1 - r_tilde)^2 + 4*r_tilde);
                else
                    kon_ij = kon0 * exp(-(dist^2)/(sigma^2));
                end

                % ----- CLAMP (CRITICAL FOR STABILITY) -----
                kon_ij = max(0, kon_ij);
                kon_ij = min(kon_ij, 1/dt);

                % ----- BINDING (ONLY IF BOTH FREE) -----
                if B(i,j) == 0 && bound(i)==0 && bound(j)==0
                    if rand < kon_ij * dt
                        B(i,j) = 1;
                        B(j,i) = 1;
                        bound(i) = 1;
                        bound(j) = 1;
                    end

                % ----- UNBINDING -----
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

    %% ----- FORCE INITIALIZATION -----
    F = zeros(N,1);

    %% ----- WLC FORCES -----
    r = diff(x);
    rabs = abs(r);
    FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);
    
    F(2:N)   = F(2:N)   - sign(r).*FWLC;
    F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

    %% ----- CROSSLINK FORCES -----
    dx = x' - x;
    Fcross = k_cross * B .* dx;
    F = F + sum(Fcross,2);

    %% ----- EXCLUDED VOLUME -----
    dx = x - x';
    FEV = cEV * dx .* exp(-aEV * dx.^2);
    FEV = sum(FEV,2);
    F = F + FEV;

    %% ----- BROWNIAN NOISE -----
    xi = noise * randn(N,1);

    %% ----- UPDATE -----
    x = x + (dt/zeta)*F + xi;

    %% ----- TETHERS -----
    x(1) = 0;
    x(N) = L;

    %% ----- STORE TRAJECTORY -----
    traj(:,t) = x;

    %% ----- STORE NETWORK SNAPSHOT -----
    if mod(t, store_interval) == 0
        B_history{snap_idx} = B;
        snap_idx = snap_idx + 1;
    end

end
toc

%% ---------------- CLUSTER ANALYSIS (PROXIMITY-BASED) ----------------
num_snaps = length(B_history);
cluster_sizes = zeros(num_snaps, N);
largest_cluster = zeros(num_snaps,1);

r_thresh = 20;   % proximity threshold (nm)

for k = 1:num_snaps
    
    x_snap = traj(:, k*store_interval);

    dx = abs(x_snap - x_snap');
    A = (dx < r_thresh) & (dx > 0);   % adjacency matrix

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

% --- 1. Crosslink matrices (actual bonds) ---
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_history{k})
    axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Crosslink Matrices (Actual Bonds)')

% --- 2. Proximity-based graphs ---
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
sgtitle('Proximity-Based Network Graphs')

% --- 3. Largest cluster over time ---
figure
plot(largest_cluster,'LineWidth',2)
xlabel('Snapshot index')
ylabel('Largest cluster size')
title('Cluster Growth (Proximity-Based)')
grid on

% --- 4. Cluster size distribution ---
all_clusters = cluster_sizes(:);
all_clusters = all_clusters(all_clusters > 0);

figure
histogram(all_clusters)
xlabel('Cluster size')
ylabel('Frequency')
title('Cluster Size Distribution')

%% ---------------- FINAL CONFIG SNAPSHOT ----------------
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

title('Final Polymer Configuration with Bonds')
ylim([-1 1])