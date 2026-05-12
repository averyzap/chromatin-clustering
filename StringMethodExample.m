% Parameters for Mueller Potential
A = [-200, -100, -170, 15];
a = [-1, -1, -6.5, 0.7];
b = [0, 0, 11, 0.6];
c = [-10, -10, -6.5, 0.7];
x0 = [1, 0, -0.5, -1];
y0 = [0, 0.5, 1.5, 1];

% String Parameters
N = 20;             % Number of points (images)
dt = 1e-4;          % Time step for steepest descent
max_iter = 5000;    % Number of iterations
phi = [zeros(N,1), linspace(1.5, -0.5, N)']; % Initial string (vertical line)

for step = 1:max_iter
    % --- Step 1: Evolution (Forward Euler) ---
    for i = 1:N
        gradV = compute_gradient(phi(i,1), phi(i,2), A, a, b, c, x0, y0);
        phi(i,:) = phi(i,:) - dt * gradV;
    end
    
    % --- Step 2: Reparameterization (Linear Interpolation) ---
    % Calculate normalized arclength (s_i)
    dist = sqrt(diff(phi(:,1)).^2 + diff(phi(:,2)).^2);
    s = [0; cumsum(dist)];
    s = s / s(end); % Normalize arclength to [0, 1]
    
    % Interpolate back to uniform arclength points (alpha_i)
    alpha = linspace(0, 1, N)';
    phi(:,1) = interp1(s, phi(:,1), alpha, 'linear');
    phi(:,2) = interp1(s, phi(:,2), alpha, 'linear');
    
    % Visualization every 500 steps
    if mod(step, 500) == 0
        plot(phi(:,1), phi(:,2), '-o'); drawnow;
    end
end

function g = compute_gradient(x, y, A, a, b, c, x0, y0)
    dx = 0; dy = 0;
    for i = 1:4
        exp_term = A(i) * exp(a(i)*(x-x0(i))^2 + b(i)*(x-x0(i))*(y-y0(i)) + c(i)*(y-y0(i))^2);
        dx = dx + exp_term * (2*a(i)*(x-x0(i)) + b(i)*(y-y0(i)));
        dy = dy + exp_term * (b(i)*(x-x0(i)) + 2*c(i)*(y-y0(i)));
    end
    g = [dx, dy];
end