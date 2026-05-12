%% TETHERED POLYMER SIMULATION (1D) WITH CROSSLINKING
clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N = 6;                 
L = 80;               
steps = 300000;
dt = 1e-3;

zeta = 2.5e-3;         
kBT  = 4.1;            

Lp = 50;
Nk = 17;
R0 = Nk*(2*Lp);        
alpha = 0.2176;        

% Excluded volume parameters
cEV = 8.305e-5;
aEV = 3.268e-5;

% Crosslink parameters
kon0  = 1.0;      % base binding rate
koff  = 0.5;      % unbinding rate
sigma = 20;       % interaction length scale (nm)
k_cross = 0.05;   % crosslink spring strength

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONDITION ----------------
x = linspace(0,L,N)';
x(1) = 0;
x(N) = L;

traj = zeros(N,steps);

% Crosslink matrix
B = zeros(N);

%% ---------------- SIMULATION ----------------
tic
for t = 1:steps

    %% ----- CROSSLINK DYNAMICS (Markov) -----
    for i = 1:N
        for j = i+1:N
            
            if abs(i-j) > 1   % skip chain neighbors
                
                dist = x(i) - x(j);
                
                % Distance-dependent binding rate
                kon_ij = kon0 * exp(-(dist^2)/(sigma^2));
                
                if B(i,j) == 0
                    % Attempt binding
                    if rand < kon_ij * dt
                        B(i,j) = 1;
                        B(j,i) = 1;
                    end
                else
                    % Attempt unbinding
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

    %% ----- CROSSLINK FORCES (FIXED) -----
    dx = x' - x;   % <-- FIXED SIGN

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

    %% ----- STORE -----
    traj(:,t) = x;
end
toc

%% ---------------- REMOVE TRANSIENT ----------------
burn_in = round(0.3*steps);
traj_ss = traj(:,burn_in:end);
T_ss = size(traj_ss,2);

%% ---------------- ANALYSIS ----------------

% 1. Mean position
mean_pos = mean(traj_ss,2);
figure
plot(1:N, mean_pos,'o-','LineWidth',2)
xlabel('Bead index')
ylabel('<x_i>')
title('Average bead positions')

% 2. Mean spacing
mean_spacing = diff(mean_pos);
figure
plot(1:N-1, mean_spacing,'o-','LineWidth',2)
xlabel('Bond index')
ylabel('<x_{i+1} - x_i>')
title('Mean bond lengths')

% 3. Variance
var_pos = var(traj_ss,0,2);
figure
plot(1:N, var_pos,'o-','LineWidth',2)
xlabel('Bead index')
ylabel('Var(x_i)')
title('Variance of bead positions')

% 4. Range explored
xmin = min(traj_ss,[],2);
xmax = max(traj_ss,[],2);
range = xmax - xmin;
figure
plot(1:N, range,'o-','LineWidth',2)
xlabel('Bead index')
ylabel('Range (nm)')
title('Spatial range per bead')

% 5. Running time average (middle bead)
running_avg = zeros(1,T_ss);
mid = ceil(N/2);
for t = 1:T_ss
   running_avg(t) = mean(traj_ss(mid,1:t));
end
figure
plot(running_avg,'LineWidth',2)
xlabel('Time')
ylabel('Running avg position')
title('Convergence of time average (middle bead)')

% 6. Histogram (middle bead)
figure
histogram(traj_ss(mid,:),50)
xlabel('Position (nm)')
title('Distribution of middle bead')

%% ---------------- OPTIONAL: TRAJECTORY PLOT ----------------
figure
plot(traj')
xlabel('Time step')
ylabel('Position (nm)')
title('Bead Position Fluctuations')