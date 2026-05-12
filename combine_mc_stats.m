%%% Load the data 

addpath('./data/');
files = dir('./data/mc_stats3_system2*.mat');

eps_range = [];
mean_times = [];

for i=1:length(files)
    fn = files(i).name;
    load(fn, 'eps','avg_time');
    eps_range(end+1) = eps;
    mean_times(end+1) = avg_time;
end

fit_x = 1./eps_range; fit_y = log(mean_times);
scatter(fit_x, fit_y);

%%% Create a fit
p1 = polyfit(fit_x(1:10), fit_y(1:10), 1);
%p1(1) = 0.015;
%p1(2) = 10.65;
p2 = polyfit(fit_x(end-8:end), fit_y(end-8:end), 1);
plot(fit_x, fit_y, 'x', fit_x, p1(1)*fit_x+p1(2),'r-')
axis([0 200 8 14]);
%legend('