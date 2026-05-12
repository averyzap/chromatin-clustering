%% 2D TETHERED POLYMER WITH WLC + EV + NOISE
clear; clc; close all;

%% ---------------- PARAMETERS ----------------

N = 25;          % number of beads
L = 400;         % wall separation (nm)

steps = 40000;
dt = 1e-3;

zeta = 2.5e-3;
kBT  = 4.1;

Lp = 50;
Nk = 17;
R0 = Nk*(2*Lp);

alpha = 0.2176;

% excluded volume
cEV = 8.305e-5;
aEV = 3.268e-5;

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL CONFIGURATION ----------------

X = zeros(N,2);

X(:,1) = linspace(0,L,N);      % x positions
X(:,2) = 20*randn(N,1);        % small random y displacement

% tether ends to walls
X(1,:) = [0 0];
X(N,:) = [L 0];

traj = zeros(N,2,steps);

%% ---------------- PLOT SETUP ----------------

figure('color','w')
hold on

xline(0,'k','LineWidth',2)
xline(L,'k','LineWidth',2)

chain = plot(X(:,1),X(:,2),'o-','LineWidth',2,'MarkerSize',6);

axis equal
xlim([-50 L+50])
ylim([-200 200])

xlabel('x (nm)')
ylabel('y (nm)')
title('2D Tethered Chromatin Chain')

drawnow

%% ---------------- SIMULATION ----------------

for t = 1:steps

    F = zeros(N,2);

    %% ---------- WLC NEIGHBOR FORCES ----------

    for i = 1:N-1

        rij = X(i+1,:) - X(i,:);
        r = norm(rij);

        if r == 0
            continue
        end

        FWLC = alpha * (-1 + 1/(1-r/R0)^2 + 4*r/R0);

        fvec = FWLC * (rij/r);

        F(i,:)   = F(i,:) + fvec;
        F(i+1,:) = F(i+1,:) - fvec;

    end

    %% ---------- EXCLUDED VOLUME ----------

    for i = 1:N
        for j = i+1:N

            rij = X(i,:) - X(j,:);
            r2 = sum(rij.^2);

            f = cEV * rij * exp(-aEV*r2);

            F(i,:) = F(i,:) + f;
            F(j,:) = F(j,:) - f;

        end
    end

    %% ---------- BROWNIAN FORCE ----------

    xi = noise * randn(N,2);

    %% ---------- UPDATE ----------

    X = X + (dt/zeta)*F + xi;

    %% ---------- ENFORCE TETHERS ----------

    X(1,:) = [0 0];
    X(N,:) = [L 0];

    traj(:,:,t) = X;

    %% ---------- ANIMATION ----------

    if mod(t,50)==0
        set(chain,'XData',X(:,1),'YData',X(:,2))
        drawnow
    end

end