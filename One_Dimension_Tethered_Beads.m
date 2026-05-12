%% TETHERED POLYMER SIMULATION (1D) WITH ANIMATION
clear; clc; close all;

%% ---------------- PARAMETERS ----------------

N = 6;                 % number of beads
L = 120;                % wall separation (nm)

steps = 30000;
dt = 1e-3;

zeta = 2.5e-3;          % drag coefficient (pN s / nm)
kBT  = 4.1;             % thermal energy

Lp = 50;
Nk = 17;
R0 = Nk*(2*Lp);         % WLC contour length

alpha = 0.2176;         % WLC force scale

% Excluded volume parameters
cEV = 8.305e-5;
aEV = 3.268e-5;

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONDITION ----------------

x = linspace(0,L,N)';   % evenly spaced start

x(1) = 0;               % left tether
x(N) = L;               % right tether

traj = zeros(N,steps);

%% ---------------- FIGURE FOR ANIMATION ----------------

figure('Color','w')
hold on

wall1 = xline(0,'k','LineWidth',2);
wall2 = xline(L,'k','LineWidth',2);

chain_plot = plot(x,1:N,'o-','LineWidth',2,'MarkerSize',6);

set(gca,'YDir','reverse')
ylim([0 N+1])
xlabel('Position (nm)')
ylabel('Bead index')
title('Tethered Polymer Dynamics')

drawnow

%% ---------------- SIMULATION ----------------

for t = 1:steps

    F = zeros(N,1);

    %% ----- WLC FORCES (nearest neighbors) -----

    r = diff(x);
    rabs = abs(r);

    FWLC = alpha * (-1 + 1./(1 - rabs/R0).^2 + 4*rabs/R0);

    % apply to beads
    F(2:N)   = F(2:N)   - sign(r).*FWLC;
    F(1:N-1) = F(1:N-1) + sign(r).*FWLC;

    %% ----- EXCLUDED VOLUME (vectorized pairwise) -----

    dx = x - x';

    FEV = cEV * dx .* exp(-aEV * dx.^2);

    FEV = sum(FEV,2);

    F = F + FEV;

    %% ----- BROWNIAN NOISE -----

    xi = noise * randn(N,1);

    %% ----- OVERDAMPED UPDATE -----

    x = x + (dt/zeta)*F + xi;

    %% ----- ENFORCE TETHERS -----

    x(1) = 0;
    x(N) = L;

    traj(:,t) = x;

    %% ----- ANIMATION -----

    if mod(t,50)==0
        set(chain_plot,'XData',x)
        drawnow
    end

end

%% ---------------- TRAJECTORY PLOT ----------------

figure
plot(traj')
xlabel('Time step')
ylabel('Position (nm)')
title('Bead Position Fluctuations')