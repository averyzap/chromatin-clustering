%% TETHERED POLYMER SIMULATION (1D)
%% FIXED N = 6, SWEEP WALL SEPARATION L
clear; clc; close all;

%% ---------------- PARAMETERS ----------------
N = 6;                         % FIXED number of beads
L_values = [150 200 300 400 600 800];   % sweep wall separation

steps = 30000;
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

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- STORAGE ----------------
mean_spacing_all = zeros(length(L_values), N-1);
var_all = zeros(length(L_values), N);
range_all = zeros(length(L_values), N);

%% ================= MAIN LOOP OVER L =================
for k = 1:length(L_values)

    L = L_values(k);

    %% INITIAL CONDITION
    x = linspace(0,L,N)';
    x(1) = 0;
    x(N) = L;

    traj = zeros(N,steps);

    %% -------- SIMULATION --------
    for t = 1:steps

        F = zeros(N,1);

        %% WLC
        r = diff(x);
        rabs = abs(r);
        FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);

        F(2:N)   = F(2:N)   - sign(r).*FWLC;
        F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

        %% EXCLUDED VOLUME
        dx = x - x';
        FEV = cEV * dx .* exp(-aEV * dx.^2);
        FEV = sum(FEV,2);

        F = F + FEV;

        %% NOISE
        xi = noise * randn(N,1);

        %% UPDATE
        x = x + (dt/zeta)*F + xi;

        %% TETHERS
        x(1) = 0;
        x(N) = L;

        traj(:,t) = x;
    end

    %% -------- REMOVE TRANSIENT --------
    burn_in = round(0.3*steps);
    traj_ss = traj(:,burn_in:end);

    %% -------- ANALYSIS --------

    % Mean position
    mean_pos = mean(traj_ss,2);

    % Mean spacing
    mean_spacing = diff(mean_pos);

    % Variance
    var_pos = var(traj_ss,0,2);

    % Range
    xmin = min(traj_ss,[],2);
    xmax = max(traj_ss,[],2);
    range = xmax - xmin;

    %% STORE
    mean_spacing_all(k,:) = mean_spacing;
    var_all(k,:) = var_pos;
    range_all(k,:) = range;

end

%% ================= PLOTTING =================

% 1. Mean spacing vs L
figure
plot(L_values, mean_spacing_all,'o-','LineWidth',2)
xlabel('Wall separation L (nm)')
ylabel('<x_{i+1} - x_i>')
title('Mean spacing vs confinement')
legend(arrayfun(@(i) sprintf('Bond %d',i),1:N-1,'UniformOutput',false))

% 2. Variance vs L
figure
plot(L_values, var_all,'o-','LineWidth',2)
xlabel('Wall separation L (nm)')
ylabel('Var(x_i)')
title('Variance vs confinement')
legend(arrayfun(@(i) sprintf('Bead %d',i),1:N,'UniformOutput',false))

% 3. Range vs L
figure
plot(L_values, range_all,'o-','LineWidth',2)
xlabel('Wall separation L (nm)')
ylabel('Range (nm)')
title('Exploration range vs confinement')
legend(arrayfun(@(i) sprintf('Bead %d',i),1:N,'UniformOutput',false))

%% ================= INTERPRETATION GUIDE =================
% Small L: compressed chain, small spacing, low variance
% Large L: stretched chain, larger spacing, higher variance (especially middle beads)
% Middle beads should always fluctuate more than tethered ends
