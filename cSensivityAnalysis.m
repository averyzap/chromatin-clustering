% Base Parameters
alpha_base = 0.2176; 
c_base = 8.305e-5;    % Base Force Strength
R0 = 1700;
a = 3.268e-5;
r_core = 25; 
khc = 0.1;           % Stiff hard-core spring

% Sensitivity range: 50% to 150% of base Force Strength (c)
variations_c = linspace(0.5, 1.5, 50); 
r_equil_c = zeros(size(variations_c));

for i = 1:length(variations_c)
    % Vary c (Excluded Volume Strength)
    c_current = c_base * variations_c(i);
    
    % Total Force function: FWLC - FEV - FHC
    F_tot = @(r) (alpha_base * (-1 + 1./(1 - r/R0).^2 + 4*r/R0)) ... % WLC [cite: 340]
                 - (c_current * r .* exp(-a * r.^2)) ...           % EV [cite: 343]
                 - (max(0, khc * (r_core - r)));                   % Hard Core [cite: 352]
             
    % Find equilibrium r*
    r_equil_c(i) = fzero(F_tot, [5, 60]); 
end

% Plotting Sensitivity for c
figure;
plot(variations_c * 100, r_equil_c, 'r', 'LineWidth', 2);
xlabel('Percentage of Base EV Strength c (%)');
ylabel('Equilibrium Separation r* (nm)');
title('Sensitivity Analysis: r* vs. EV Strength (c)');
grid on;