function [mean_times, all_times] = escape_statistics(system, x, eps_range, trials)
x = find_minima(system, x);
mean_times = zeros(size(eps_range));
all_times = cell(size(eps_range));
x2r = reshape(x, system.d, []);

dtmax = 10000000;

if nargin < 4
    trials = 5000;
end

parfor (e_c = 1:length(eps_range), 4)
    eps = eps_range(e_c);
    %disp(eps);
    times = zeros(1, trials);
    %h = waitbar(0, 'Simulating');
    for trial=1:trials
        %disp(trial)
        [~, ~, time] = simulate_system(system, x2r, 0.001, dtmax, eps,1);
        times(trial) = time;
        %waitbar(trial/trials, h);
        if mod(trial, 20) ==0
            disp(trial);
        end
    end
    %close(h);
    avg_time = sum(times)/sum(times~=(dtmax));
    %disp(avg_time)
    mean_times(e_c) = avg_time;
    all_times{e_c} = times;
    %hist(times) 
end
end