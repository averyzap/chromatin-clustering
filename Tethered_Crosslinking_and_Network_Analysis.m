%% TETHERED POLYMER SIMULATION (1D) WITH CROSSLINKING + NETWORK ANALYSIS
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
k_cross = 0.01;   % safer value

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONDITION ----------------
x = linspace(0,L,N)';
x(1) = 0;
x(N) = L;

traj = zeros(N,steps);

% Crosslink matrix
B = zeros(N);

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
                kon_ij = kon0 * exp(-(dist^2)/(sigma^2));
                
                if B(i,j) == 0
                    if rand < kon_ij * dt
                        B(i,j) = 1;
                        B(j,i) = 1;
                    end
                else
                    if rand < koff * dt
                        B(i,j) = 0;
                        B(j,i) = 0;
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

    %% ----- CROSSLINK FORCES (CORRECT SIGN) -----
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

%% ---------------- CLUSTER ANALYSIS ----------------
num_snaps = length(B_history);
cluster_sizes = zeros(num_snaps, N);
largest_cluster = zeros(num_snaps,1);

for k = 1:num_snaps
    
    Bsnap = B_history{k};
    G = graph(Bsnap);
    
    bins = conncomp(G);
    
    sizes = zeros(1,max(bins));
    for c = 1:max(bins)
        sizes(c) = sum(bins == c);
    end
    
    cluster_sizes(k,1:length(sizes)) = sizes;
    largest_cluster(k) = max(sizes);
end

%% ---------------- VISUALIZATION ----------------

% --- 1. Crosslink matrices ---
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    imagesc(B_history{k})
    axis square
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Crosslink Matrices Over Time')

% --- 2. Graph visualization ---
figure
for k = 1:num_snaps
    subplot(ceil(sqrt(num_snaps)), ceil(sqrt(num_snaps)), k)
    G = graph(B_history{k});
    plot(G,'Layout','force')
    title(['t = ', num2str(k*store_interval)])
end
sgtitle('Crosslink Network Graphs')

% --- 3. Largest cluster over time ---
figure
plot(largest_cluster,'LineWidth',2)
xlabel('Snapshot index')
ylabel('Largest cluster size')
title('Cluster Growth Over Time')
grid on

% --- 4. Cluster size distribution ---
all_clusters = cluster_sizes(:);
all_clusters = all_clusters(all_clusters > 0);

figure
histogram(all_clusters)
xlabel('Cluster size')
ylabel('Frequency')
title('Cluster Size Distribution')

%% ---------------- OPTIONAL: FINAL CONFIG SNAPSHOT ----------------
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

title('Final Polymer Configuration with Crosslinks')
ylim([-1 1])