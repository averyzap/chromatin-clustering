quasistrings = {};
h = waitbar(0, 'Simulating strings');
for i=1:length(split_ic)
    ic = split_ic{i};
    out = string_v2(system2, ic, 1e-3, 100, 'maxiters',500);
    
    % figure out if it is monotonic
    quasistrings{i} = out;
    
    waitbar(i/length(strings2), h);
end
close(h)