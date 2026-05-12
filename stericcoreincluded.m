%% ============================================================
%  WLC + Excluded Volume + Hard-Core Repulsion
%  Corrected Model (No Collapse)
% ============================================================

clear; clc; close all;

%% -------------------------------
% Physical constants
% -------------------------------
kBT  = 4.1;                 % pN*nm
zeta = 2.5e-3;              % pN*s/nm

%% -------------------------------
% WLC parameters
% -------------------------------
R0    = 1700;               % nm
alpha = 0.2176;             % pN

%% -------------------------------
% Gaussian EV parameters (Ben-style)
% -------------------------------
cEV = 8.305e-5;             % pN / nm
aEV = 3.268e-5;             % nm^-2

%% -------------------------------
% Hard-core steric repulsion
% -------------------------------
r_core = 25;                % nm (steric bead radius)
k_hc   = 5e-3;              % pN / nm  (>> WLC stiffness)

%% ============================================================
%  SECTION 1: FORCE CURVES
% ============================================================

r = linspace(0,60,3000);

FWLC = WLC_force(r,alpha,R0);
FEV  = EV_force(r,cEV,aEV);
FHC  = HC_force(r,r_core,k_hc);

FREP = FEV + FHC;           % total repulsion
G    = FWLC - FREP;

% Find equilibrium
idx = find(diff(sign(G)),1);
r_star = interp1(G(idx:idx+1), r(idx:idx+1), 0);

figure;
plot(r,FWLC,'LineWidth',2); hold on;
plot(r,FREP,'LineWidth',2);
xline(r_star,'k--','LineWidth',1.5);
xlabel('Separation r (nm)');
ylabel('Force (pN)');
legend('WLC attraction','Total repulsion (EV + hard core)','Equilibrium');
title('Force Balance with Hard-Core Repulsion');
grid on;

fprintf('Equilibrium separation r* ≈ %.2f nm\n',r_star);

%% ============================================================
%  SECTION 2: 1D OVERDAMPED DYNAMICS
% ============================================================

dt   = 1e-3;        % s
Tend = 10;          % s
Nt   = round(Tend/dt);
t    = linspace(0,Tend,Nt);

r0_list = [10, 50];
colors = lines(2);

figure; hold on;

for k = 1:length(r0_list)
    r_traj = zeros(1,Nt);
    r_traj(1) = r0_list(k);

    for n = 1:Nt-1
        Ftot = WLC_force(r_traj(n),alpha,R0) ...
             - EV_force(r_traj(n),cEV,aEV) ...
             - HC_force(r_traj(n),r_core,k_hc);

        r_traj(n+1) = r_traj(n) - dt*Ftot/zeta;
        r_traj(n+1) = max(r_traj(n+1),0);
    end

    plot(t,r_traj,'LineWidth',2,...
        'DisplayName',['r(0) = ' num2str(r0_list(k)) ' nm']);
end

yline(r_star,'k--','Equilibrium');
xlabel('Time (s)');
ylabel('Separation r(t) (nm)');
title('Stable Overdamped Dynamics with Steric Core');
legend;
grid on;

%% ============================================================
%  LOCAL FORCE FUNCTIONS
% ============================================================

function F = WLC_force(r,alpha,R0)
    x = r./R0;
    F = alpha.*(-1 + 1./(1-x).^2 + 4*x);
end

function F = EV_force(r,c,a)
    F = c.*r.*exp(-a.*r.^2);
end

function F = HC_force(r,r_core,k)
    F = zeros(size(r));
    mask = r < r_core;
    F(mask) = k*(r_core - r(mask));
end
