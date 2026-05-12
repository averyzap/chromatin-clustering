function outstruct = string_descend(system, x, dt, n1, varargin)
    % Document function

    if ndims(x) == 3
        x = reshape(x, [], size(x, 3));
    end

    p = inputParser;
    addOptional(p, 'checkpoint_at', -1);
    addOptional(p, 'plot_during', 0);
    addOptional(p, 'make_movie', 0);
    addOptional(p, 'nn', 1);
    addOptional(p, 'tol', 1e-6);
    addOptional(p, 'climbing', 0);
    addOptional(p, 'fix_beads', 1);
    addOptional(p, 'maxiters', 10000);
    addOptional(p, 'deterministic', 0);
    addOptional(p, 'dW', zeros(size(x,1),n1));
    addOptional(p, 'verbose', 0);
    addOptional(p, 'implicit', 1);

    parse(p, varargin{:});

    % set up input options
    checkpoint_at = p.Results.checkpoint_at;
    plot_during= p.Results.plot_during;
    make_movie = p.Results.make_movie;
    nn = p.Results.nn;   % when to plot
    diff_tol1 = p.Results.tol;
    climbing = p.Results.climbing;
    fix_beads = p.Results.fix_beads;
    nstepmax = p.Results.maxiters;
    deterministic = p.Results.deterministic;
    dW = p.Results.dW;
    verbose = p.Results.verbose;
    implicit = p.Results.implicit;

    if fix_beads
        % rotate system to maintain symmetry fixing
        for j=1:size(x, 2)
            xv = x(5:6, j);
            r = norm(xv);
            c = xv(2)/r;
            s = xv(1)/r;
            R = [c -s; s c];
            x(1:2, j) = R*x(1:2, j);
            x(3:4, j) = R*x(3:4, j);
            x(5:6, j) = R*x(5:6, j);
        end
    end

    % Set up for file name in checkpoints
    ts = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss'));

    d = system.d;
    p = size(x,1);
    g1 = linspace(0,1,n1);

    % redistribute images along string at constant arclength
    dx = diff(x,1,2);

    dd = sum(dx.^2);
    dd = sqrt([0 dd]);

    ll = cumsum(dd);
    ll = ll/ll(end);
    x = interp1(ll',x',g1', 'makima')';

    [W, S] = compute_W(system, x);

    xOld = x;

    x_arr = {};

    if plot_during || make_movie
        bounds = repmat([-1 1],[1 system.d]);
    end

    if make_movie
        [~, msg] = mkdir('movies');
        v = VideoWriter(['movies/' ts '.mp4'], 'MPEG-4');
        open(v);
        frames = [];
    end


    exitflag = 0;   % changes to 1 if convergence reached
    Sincflag = deterministic;    % don't give this warning for deterministic

    diffX = zeros(size(x));

    W_arr = zeros(1, nstepmax);
    err_arr = zeros(1, nstepmax);

    if climbing
        end_val = n1;
    else
        end_val = n1-1;
    end

    if ~climbing
        dW(:, 1) = 0;
        dW(:, end) = 0;
    end

    S = nan;
    for nstep = 1:nstepmax
        x0 = x;

        % Compute dphi using the interpolation scheme
        dphi = x(:, 3:end)-x(:, 1:(end-2));
        dphi = [dphi(:, 1) dphi dphi(:, end)];

        if implicit
            % backward euler implicit loop
            xnew = x;
            sub_err = [0 ones(1,n1-2) 0];
            sub_iter = 0;
            dW_pred = dW;
            xnew = x;
            be_tol = 1e-15;
            while 1
                dphi = xnew(:, 3:end)-xnew(:, 1:(end-2));
                dphi = dphi ./ vecnorm(dphi);
                dphi = [dphi(:, 1) dphi dphi(:, end)];
                dWnew = dW_sub(system, xnew, dphi, dW, deterministic, find(sub_err>be_tol));
                dW = dWnew;
                if sub_iter > 0
                    sub_diff = (x-dt*dWnew) - xnew;
                    sub_err = vecnorm(sub_diff);
                    if max(sub_err) < be_tol
                        break
                    end
                end
                xnew = x-dt*dW;


                sub_iter = sub_iter + 1;
                if sub_iter > 20
                    asdf = 1;
                end
            end
        else
            dWnew = dW_sub(system, x, dphi, dW, deterministic);
            xnew = x - dt*dWnew;
        end

        % climbing doesn't use the backward euler
        if climbing
            dW(:, end) =  compute_dW(system, x(:,end),  dphi(:,end),  dW(:,end));
        end

        Sold = S;
        [Wmax, S, W] = getW(x,dW);

        if climbing
            % change final image to have the climbing direction
            change = dW(:, end) - 1.5*dphi(:,end)*dot(dW(:,end), dphi(:, end))/(norm(dphi(:,end))^2);
            change = -dt*change;
            xnew(:, end) = x(:, end) + change;
        end

        x = xnew;
        dW = dWnew;


        if fix_beads
            % rotate system to maintain symmetry fixing
            for j=1:end_val
                xv = x(5:6, j);
                r = norm(xv);
                c = xv(2)/r;
                s = xv(1)/r;
                R = [c -s; s c];
                x(1:2, j) = R*x(1:2, j);
                x(3:4, j) = R*x(3:4, j);
                x(5:6, j) = R*x(5:6, j);
            end
        end


        % redistribute along string
        dx = x(:,2:n1) - x(:,1:n1-1);

        dd = sum(dx.^2);
        dd = sqrt([0 dd]);

        ll = cumsum(dd);
        string_len = ll(end);
        ll = ll/ll(n1);
        xnew = interp1(ll',x',g1', 'linear')';
        x=xnew;


        % --- exit condition
        diffold = diffX;
        diffX = x-x0;
        diffX_arr = diffX(:);
        diffFromStart = x-xOld;


        errX = norm(diffX(:), Inf);
        errX2 = norm(diffX(:), 1);
        errFromStart = norm(diffFromStart(:), Inf);
        errEnd = norm(diffX(:, end), Inf);
        err_align = dot(diffold(:), diffX(:))/(norm(diffold(:))*norm(diffX(:)));

        if verbose && mod(nstep,100)==1
            fprintf("Iteration: %d  W max:%.3e  Error:%.3e  S: %e  Alignment:%f Length: %e\n", nstep,Wmax, errX, S, err_align, string_len);
        end
        W_arr(nstep) = Wmax;
        err_arr(nstep) = errX;

        x_arr{numel(x_arr)+1} = x;

        if nstep == 500
            xOld = x;
        end

        if (make_movie || plot_during) && mod(nstep, nn) == 0
            subplot(121);
            plot(g1,x');
            subplot(122);
            dx = x(:, 2:end)-x(:, 1:(end-1));
            dx = [dx dx(:, end)];

            plot(g1,dx');
            F = getframe;   % also forces update of plot
            if make_movie
                writeVideo(v, F);
            end
        end

         if checkpoint_at > 0 && mod(nstep, checkpoint_at)==0
            outstruct = [];
            outstruct.x = x;
            outstruct.dW = dW;
            outstruct.exitflag = exitflag;
            outstruct.W_arr = W_arr;
            outstruct.err_arr = err_arr;
            outstruct.W = W;

            % if you save this the checkpoint files are extremely large
            %outstruct.x_arr = x_arr;

            out = outstruct;
            save(['checkpoints/' ts '_' num2str(nstep) '.mat'], 'out');
        end
        if errX < diff_tol1*dt || (1-Wmax/S) < 1e-4
            exitflag = 1;
            W_arr = W_arr(1:nstep);
            err_arr = err_arr(1:nstep);
            break
        end
    end

    if make_movie
        close(v);
        close all;
    end

    outstruct = [];
    outstruct.x = x;
    outstruct.dW = dW;
    outstruct.exitflag = exitflag;
    outstruct.W_arr = W_arr;
    outstruct.err_arr = err_arr;
    outstruct.x_arr = x_arr;
    outstruct.W = W;
end

function [Wmax, S, W] = getW(x, dW)
    dx = (x(:,2:end) - x(:,1:(end-1)));
    % fill in dx at endpoints
    dx = [dx x(:,end)-x(:,(end-1))];
    %W = cumsum(sum(dW.*dx));
    dWp = dW;
    W = cumsum(sum(0.5*(dWp(:, 2:end).*dx(:, 2:end) + dWp(:, 1:(end-1)).*dx(:,1:(end-1)))));
    W2 = vecnorm(dW).*vecnorm(dx);
    W2(dot(dW, dx)<=0) = 0;
    W2 = cumsum(W2);
    W = [0 W];
    [Wmax] = max(W);
    S = max(W2);
end

function [dW, conv_flag, dU, d2W] = dW_sub(system, x, dphi, guess, deterministic, images)
    dW = guess;
    dU = zeros(size(guess));
    if nargout == 4
        d2W = zeros(size(guess));
    end
    n1 = size(guess, 2);
    if deterministic
        conv_flag = 0;
    end
    if nargin < 6
        images = 2:(n1-1);
    end
    for j=images
        guess = dW(:, j);
        
        % figure out which direction we should compute the quasipotential
        % in
        dg = deterministic_force(system, x(:, j));
        dU(:, j) = dg;
        if deterministic
            dW(:, j) = -dg;
            continue;
        end
        dp = dot(dg, dphi(:, j));

        [tmp, conv_flag(j-1)] =  compute_dW(system, x(:,j),  -sign(dp)*dphi(:,j),  guess);

        d2w_flag = 0; 
        if nargout == 4
            eps = 1e-7;
            [d2W(:, j), d2w_flag] = compute_dW(system, x(:, j) + eps*tmp, -sign(dp)*dphi(:,j), tmp);
             d2W(:, j) = (d2W(:, j) - tmp)/eps;
        end
        
        if conv_flag(j-1)~=0 || d2w_flag ~= 0
            fprintf("Failed image %d\n",j);
            dW(:, j) = 0;
        else
            dW(:,j) = tmp;
        end
    end
end