function plot2dbeads(pos, ax)
    if ndims(pos) == 2 && size(pos,1)>2
        pos = reshape(pos, 2, [], size(pos,2));
    end
    if nargin < 2
        ax = gca;
    end
    [d, p, nt] = size(pos);
    
%     if nargin < 2
%         figure;
%     %end
    if nt > 1
        % plotting trajectory
        cla;
        hold on;
        %cm = lines(p);
        cm = hsv(p);
        cm(1, :) = [1 0 1];
        %cm([2 3],:) = cm([3 2],:);
        cmap = colormap(cm);
        if nt > 5000
            lw = 1;
        else
            lw = 2.5;
        end
        for k=1:p
            plot(squeeze(pos(1,k,:)), squeeze(pos(2,k,:)), 'LineWidth', lw, 'Color', cmap(k,:));
        end
        hold off;
        axis([-1 1 -1 1])
    else
        
        % just plotting a single set of positions
        c=linspace(1, 5, p)';
        %ccm = [colormap(hsv(p)); colormap(lines(p))]
        %colormap(ccm);
        
        colormap(lines)
        plot(pos(1,[1 2],1), pos(2, [1 2],1), pos(1,[1 3],1), pos(2, [1 3],1), pos(1,[2 3],1), pos(2, [2 3],1), 'LineWidth', 3);
        
        hold on;
        c = hsv(p);
        c(1, :) = [1 0 1];
        pv = scatter(ax, squeeze(pos(1,:,1)), squeeze(pos(2,:,1)), 125, colormap(c), 'filled');
        
        hold off;
        axis([-1 1 -1 1])
    end
end