load('systems.mat')
regen = 1;
posfile = 'data/network_pos.mat';
if regen || ~isfile(posfile)
    % bound state
    x(:,1) = find_minima(system2, [0 -1 0 -1 0 1]');
    % small line
    x(:,2) = find_minima(system2, [0 -0.5 0 0 0 0.5]');
    % small triangle
    x(:,3) = find_minima(system2, [0.2 0 -0.2 0 0 0.2]');
    save(posfile,'x');
else
    load(posfile);
end

string_ic = {};

for i=1:2
    for j=(i+1):3
        string_ic{end+1} = [x(:,i) x(:,j)];
    end
end

strings = {};
h = waitbar(0, 'Simulating strings');
for i=1:length(string_ic)
    ic = string_ic{i};
    out = string_descend(system2, ic, 1e-3, 20, 'maxiters',2000,'deterministic',1);
    
    % figure out if it is monotonic
    strings{i} = out;
    
    waitbar(i/length(string_ic), h);
end
close(h)