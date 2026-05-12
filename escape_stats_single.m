function [avg_time] = escape_stats_single(sys_name, eps, trials, eps_seed)
    if nargin < 3
        trials = 1000;
    end
    if nargin < 4
        eps_seed = 0;
    end
    % for reproducibility, seed based on epsilon
    if eps_seed~=0
        rng(round(eps_seed/eps));
    end
    S = load('systems.mat');
    system = S.(sys_name);
    x = [0 -1 0 -1 0 1];
    x = find_minima(system, x);
    x2r = reshape(x, system.d, []);

    dtmax = 10000000;

    times = zeros(1, trials);
    for trial=1:trials
        [~,~,time] = simulate_system(system, x2r, 0.001, dtmax, eps, 1);
        times(trial) = time;
        if mod(trial, 20) ==0
            disp(trial);
        end
    end
    avg_time = sum(times)/sum(times~=(dtmax));
    if nargout == 0
        save(['data/mc_stats4_' sys_name '_' num2str(eps) '.mat'], 'eps', 'avg_time','times');
    end
end
