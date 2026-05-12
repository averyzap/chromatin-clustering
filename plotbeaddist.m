% plot all pairwise bead distances as a function of time
function plotbeaddist(pos, dt)
    [d, p, nt] = size(pos);
     
    stride = 1;
    
    if nargin < 2
        kv = 1:stride:nt;
        vmax = nt;
    else
        kv = linspace(0, dt*nt, nt);
        vmax = dt*nt;
    end
    dists = zeros(p*(p-1)/2, length(kv));
    
    for c=1:length(kv)
        dists(:, c) = pdist(pos(:, :, c)');
    end
    %hold on;
    plot(kv, dists')
    if p==3
        legend('1-2', '1-3', '2-3');
    end
    
    axis([0 vmax 0 2.4]);
    xlabel('Time')
    ylabel('Distance');

   % hold off;
end