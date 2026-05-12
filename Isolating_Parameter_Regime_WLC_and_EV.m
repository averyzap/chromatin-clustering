%% ============================================================
%  WLC + Excluded Volume Parameter Regime Analysis
%  Goal: Identify parameter regime where EV force alone
%  stabilizes the WLC bond (no steric core)
%
%  Author: Avery Zapata
% =============================================================

clear; close all; clc;

%% =============================================================
% Physical parameters
% =============================================================

kBT = 4.1;                     % pN*nm
zeta = 2.5e-3;                 % pN*s/nm

alpha = 0.2176;                % WLC force scale
R0 = 1700;                     % nm

% baseline EV parameters
c_base = 8.305e-5;             % pN/nm
a_base = 3.268e-5;             % nm^-2

%% =============================================================
% Simulation grid
% =============================================================

r = linspace(0,100,2000);      % distance grid

%% =============================================================
% FORCE DEFINITIONS
% =============================================================

FWLC = @(r,alpha) alpha .* (-1 + 1./(1 - r./R0).^2 + 4*r./R0);

FEV  = @(r,c,a) c .* r .* exp(-a .* r.^2);

%% =============================================================
% Plot force curves (diagnostic)
% =============================================================

figure
plot(r, FWLC(r,alpha),'LineWidth',2)
hold on
plot(r, FEV(r,c_base,a_base),'LineWidth',2)

xlabel('Separation r (nm)')
ylabel('Force (pN)')
legend('WLC','Excluded Volume')
title('Force Comparison')
grid on

%% =============================================================
% TOTAL FORCE FUNCTION
% =============================================================

Ftot = @(r,c,a) FWLC(r,alpha) - FEV(r,c,a);

%% =============================================================
% FUNCTION TO FIND EQUILIBRIUM SEPARATION
% =============================================================

find_equilibrium = @(c,a) equilibrium_root(Ftot,c,a);

%% =============================================================
% BIFURCATION ANALYSIS
% Sweep EV strength c
% =============================================================

c_vals = logspace(-7,-2,150);
r_star = nan(size(c_vals));

for i = 1:length(c_vals)

    c = c_vals(i);

    r_guess = 30;

    try
        r_star(i) = fzero(@(x) Ftot(x,c,a_base), r_guess);
    catch
        r_star(i) = NaN;
    end

end

figure
semilogx(c_vals, r_star,'LineWidth',2)
xlabel('EV Strength c')
ylabel('Equilibrium Separation r* (nm)')
title('Bifurcation Diagram: Collapse vs Stabilization')
grid on

%% =============================================================
% 2D PARAMETER SWEEP
% (c vs a stability map)
% =============================================================

c_grid = logspace(-7,-2,80);
a_grid = logspace(-6,-4,80);

stable_map = zeros(length(c_grid), length(a_grid));

for i = 1:length(c_grid)

    for j = 1:length(a_grid)

        c = c_grid(i);
        a = a_grid(j);

        try
            root = fzero(@(x) Ftot(x,c,a),30);

            if root > 1
                stable_map(i,j) = 1;
            end

        catch
            stable_map(i,j) = 0;
        end

    end

end

figure
imagesc(log10(a_grid), log10(c_grid), stable_map)
set(gca,'YDir','normal')

xlabel('log10(a)')
ylabel('log10(c)')
title('Stability Map (1 = Stable Nonzero Equilibrium)')
colorbar

%% =============================================================
% SENSITIVITY ANALYSIS
% =============================================================

param_variation = linspace(-0.2,0.2,60);

r_sensitivity = zeros(size(param_variation));

for i = 1:length(param_variation)

    c = c_base * (1 + param_variation(i));

    try
        r_sensitivity(i) = fzero(@(x) Ftot(x,c,a_base),30);
    catch
        r_sensitivity(i) = NaN;
    end

end

figure
plot(param_variation*100, r_sensitivity,'LineWidth',2)

xlabel('Percent Change in c (%)')
ylabel('Equilibrium r* (nm)')
title('Sensitivity Analysis of EV Strength')
grid on

%% =============================================================
% TIME DYNAMICS SIMULATION
% =============================================================

dt = 0.01;
T = 200;

steps = round(T/dt);

r_sim = zeros(steps,1);
r_sim(1) = 50;

for i = 2:steps

    r_current = r_sim(i-1);

    force = Ftot(r_current,c_base,a_base);

    rdot = -force / zeta;

    r_sim(i) = r_current + dt * rdot;

end

time = (0:steps-1)*dt;

figure
plot(time, r_sim,'LineWidth',2)
xlabel('Time')
ylabel('Separation r(t)')
title('Dynamic Evolution of Bead Separation')
grid on

%% =============================================================
% MULTIPLE INITIAL CONDITIONS TEST
% =============================================================

initial_conditions = [10 20 50 80];

figure
hold on

for ic = 1:length(initial_conditions)

    r_sim = zeros(steps,1);
    r_sim(1) = initial_conditions(ic);

    for i = 2:steps

        r_current = r_sim(i-1);

        force = Ftot(r_current,c_base,a_base);

        rdot = -force / zeta;

        r_sim(i) = r_current + dt * rdot;

    end

    plot(time,r_sim,'LineWidth',2)

end

xlabel('Time')
ylabel('r(t)')
title('Convergence from Multiple Initial Conditions')
legend('10 nm','20 nm','50 nm','80 nm')
grid on

%% =============================================================
% HELPER FUNCTION
% =============================================================

function rstar = equilibrium_root(Ftot,c,a)

    try
        rstar = fzero(@(x) Ftot(x,c,a),30);

        if rstar < 1
            rstar = NaN;
        end

    catch
        rstar = NaN;
    end

end