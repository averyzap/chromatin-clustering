% Base Parameters
alpha_base = 0.2176; 
c_base = 8.305e-5;
R0 = 1700;
a = 3.268e-5;
r_core = 25; % From your 1/20 notes
khc = 0.1;   % High stiffness relative to kWLC

% Sensitivity range: 80% to 120% of base value
variations = linspace(0.8, 1.2, 50); 
r_equil = zeros(size(variations));

for i = 1:length(variations)
    % Vary Alpha (WLC Strength)
    alpha = alpha_base * variations(i);
    
    % Define the Total Force function
    F_tot = @(r) (alpha * (-1 + 1./(1 - r/R0).^2 + 4*r/R0)) ... % WLC
                 - (c_base * r .* exp(-a * r.^2)) ...          % EV
                 - (max(0, khc * (r_core - r)));               % Hard Core
             
    % Find where Force = 0 (equilibrium)
    r_equil(i) = fzero(F_tot, [5, 60]); 
end

% Plotting Sensitivity
plot(variations * 100, r_equil, 'LineWidth', 2);
xlabel('Percentage of Base Alpha (%)');
ylabel('Equilibrium Separation r* (nm)');
title('Sensitivity Analysis: r* vs. WLC Strength (\alpha)');
grid on;