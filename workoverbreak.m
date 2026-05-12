%% ============================================================
%  WLC + Excluded Volume Validation Suite
%  Avery Zapata – End-of-Break Deliverables
% =============================================================

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
% Excluded Volume parameters
% -------------------------------
cEV = 8.305e-5;             % pN / nm
aEV = 3.268e-5;             % nm^-2

%% ============================================================
%  SECTION 1: FORCE CURVES
% ============================================================

r = linspace(0,60,2000);    % nm

FWLC = WLC_force(r,alpha,R0);
FEV  = EV_force(r,cEV,aEV);

% Locate equilibrium
G = FWLC - FEV;
idx = find(diff(sign(G)),1);
r_star = interp1(G(idx:idx+1), r(idx:idx+1), 0);

figure;
plot(r,FWLC,'LineWidth',2); hold on;
plot(r,FEV,'LineWidth',2);
xline(r_star,'k--','LineWidth',1.5);
xlabel('Separation r (nm)');
ylabel('Force (pN)');
legend('WLC attraction','EV repulsion','Equilibrium');
title('Force Balance: WLC vs Excluded Volume');
grid on;

fprintf('Equilibrium separation r* ≈ %.2f nm\n',r_star);

%% ============================================================
%  SECTION 2: 1D OVERDAMPED DYNAMICS
% ============================================================

dt   = 1e-3;        % s
Tend = 10;          % s
Nt   = round(Tend/dt);
t    = linspace(0,Tend,Nt);

r0_list = [10, 50]; % initial conditions
colors = lines(2);

figure; hold on;

for k = 1:length(r0_list)
    r_traj = zeros(1,Nt);
    r_traj(1) = r0_list(k);

    for n = 1:Nt-1
        Ftot = WLC_force(r_traj(n),alpha,R0) ...
             - EV_force(r_traj(n),cEV,aEV);
        r_traj(n+1) = r_traj(n) - dt*Ftot/zeta;
        r_traj(n+1) = max(r_traj(n+1),0); % enforce r >= 0
    end

    plot(t,r_traj,'Color',colors(k,:), ...
         'LineWidth',2, ...
         'DisplayName',['r(0) = ' num2str(r0_list(k)) ' nm']);
end

yline(r_star,'k--','Equilibrium');
xlabel('Time (s)');
ylabel('Separation r(t) (nm)');
title('1D Overdamped Dynamics');
legend;
grid on;

%% ============================================================
%  SECTION 3: BIFURCATION / PARAMETER SWEEP
% ============================================================

Sk_vals = linspace(20,120,50);     % nm
nu_vals = logspace(2,7,60);

Rstar = nan(length(Sk_vals),length(nu_vals));

for i = 1:length(Sk_vals)
    Sk = Sk_vals(i);

    beta = 3/(2*Sk^2);
    A    = (3/(2*pi*Sk^2))^(3/2);
    C    = kBT*A*beta;

    for j = 1:length(nu_vals)
        nu = nu_vals(j);

        FEV_gauss = @(r) nu*C*r.*exp(-beta*r.^2);
        FWLC_fun  = @(r) WLC_force(r,alpha,R0);

        Gfun = @(r) FWLC_fun(r) - FEV_gauss(r);

        try
            r_sol = fzero(Gfun,[1 80]);
            if r_sol > 0 && r_sol < 100
                Rstar(i,j) = r_sol;
            end
        catch
            % no equilibrium
        end
    end
end

figure;
imagesc(log10(nu_vals),Sk_vals,Rstar);
set(gca,'YDir','normal');
colorbar;
xlabel('log_{10}(\nu)');
ylabel('S_k (nm)');
title('Equilibrium Separation r^* (nm)');
colormap turbo;

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function F = WLC_force(r,alpha,R0)
    x = r./R0;
    F = alpha.*(-1 + 1./(1-x).^2 + 4*x);
end

function F = EV_force(r,c,a)
    F = c.*r.*exp(-a.*r.^2);
end
