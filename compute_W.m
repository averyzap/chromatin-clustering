function [W, S, dW, U] = compute_W(system, x, interpolate, guess)
    if nargin < 4
        guess = [];
    end
    if nargin < 3
        interpolate = 0;
        n1 = size(x, 2);
    else
        n1 = 100;
    end
    if ndims(x) == 3
        x = reshape(x, [], size(x, 3));
    end
    p = size(x, 1);
    
    if interpolate
        g1 = linspace(0,1,n1);
        % redistribute along string
        dx = diff(x,1,2);

        dd = sum(dx.^2);
        dd = sqrt([0 dd]);

        ll = cumsum(dd);
        ll = ll/ll(end);
        x = interp1(ll',x',g1', 'linear')';
    end
    
    %dx = [x(:, 2:end) - x(:, 1:(end-1)) x(:, end)-x(:, end-1)];
    dphi = (x(:, 3:end)-x(:, 1:(end-2)))/2;
    dx = [dphi(:, 1) dphi dphi(:, end)];
    
    dW = zeros(p, n1);
    dW2 = dW;
    conv_flag = zeros(n1,1);
    r_arr = zeros(4, n1);
    for j = 2:(n1-1)
    %parfor (j = 2:(n1-1), 0)
        [S, r_arr(:, j)] = switching(system, x(:, j));
        
        
        V = drift_mat(system, x(:, j));
        dg = deterministic_force(system, x(:, j));
        dp = dot(dg, dx(:, j-1));
        % figure out which direction we should compute the quasipotential
        % in
        if ~isempty(guess)
            subguess = guess(:, j);
        else
            subguess = dg/2;
        end
        
        [dW(:, j), conv_flag(j)] = compute_dW(system, x(:,j), -sign(dp)*dx(:,j), subguess);
        
    end
    %W = cumsum(sum(dW.*dx));
    W = cumsum(sum(0.5*(dW(:, 2:end).*dx(:, 2:end) + dW(:, 1:(end-1)).*dx(:,1:(end-1)))));
    W = [0 W];
    
    W2 = cumsum(vecnorm(dW).*vecnorm(dx));
    S = W2(end);
    
    if nargout >= 4
        U = deterministic_U(system, x);
    end
end