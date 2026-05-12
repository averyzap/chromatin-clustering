clc; clear; close all;

%% --- Define transition rate matrix ---
% K(i,j) = rate from state i -> j
K = [ 0    1.0   0.5;
      0.2  0     0.8;
      0.3  0.4   0 ];

num_states = size(K,1);

%% --- Simulation parameters ---
T_max = 50;        % total simulation time
t = 0;             % current time
state = 1;         % initial state

time_trace = [0];
state_trace = [state];

%% --- Run CTMC simulation ---
while t < T_max
    
    % Total rate out of current state
    lambda = sum(K(state, :));
    
    if lambda == 0
        break; % absorbing state
    end
    
    % Sample waiting time (exponential)
    tau = exprnd(1/lambda);
    
    % Advance time
    t = t + tau;
    
    % Choose next state
    probs = K(state, :) / lambda;
    next_state = randsample(1:num_states, 1, true, probs);
    
    % Store
    time_trace(end+1) = t;
    state_trace(end+1) = next_state;
    
    % Update state
    state = next_state;
end

%% --- Plot state vs time (step plot) ---
figure;
stairs(time_trace, state_trace, 'LineWidth', 2);
xlabel('Time');
ylabel('State');
title('Continuous-Time Markov Chain Trajectory');
grid on;

%% --- Optional: Histogram of state occupancy ---
dt = 0.1;
t_grid = 0:dt:T_max;
state_interp = interp1(time_trace, state_trace, t_grid, 'previous');

figure;
histogram(state_interp, 'Normalization', 'probability');
xlabel('State');
ylabel('Probability');
title('State Occupancy Distribution');