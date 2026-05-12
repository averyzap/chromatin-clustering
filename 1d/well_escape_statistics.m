function [out] = well_escape_statistics(system, eps_range, dtmax, trials)
    %eps_range = linspace(0.15, 0.2, 6);
    %dtmax = 100000;
    %trials = 5;
    mean_times = zeros(size(eps_range));
    all_times = cell(size(eps_range));

    
    for e_c = 1:length(eps_range)
        eps = eps_range(e_c);
        times = zeros(1, trials);
        h = waitbar(0, 'Simulating'); 
        for trial=1:trials
            xout = simulate_switching_well(system, 0, 1e-3, dtmax, eps,1);
            times(trial) = length(xout);
            waitbar(trial/trials, h);
        end
        close(h);
        avg_time = sum(times)/sum(times~=(dtmax+1));
        %disp(avg_time)
        mean_times(e_c) = avg_time;
        all_times{e_c} = times;
        %hist(times)
    end
    out = [];
    out.system = system;
    out.mean_times = mean_times;
    out.all_times = all_times;
    out.eps_range = eps_range;
end