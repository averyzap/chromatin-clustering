function [W, H_arr, r, dW] = compute_W(system, x, guess)
    if nargin < 3
        guess = [];
    end
    if ndims(x) == 3
        x = reshape(x, [], size(x, 3));
    end
    p = size(x, 1);
    n1 = size(x, 2);
    n = system.n;
    
    dW_fun = @compute_dW_1d;
    
    % number of images along the string (try from  n1 = 3 up to n1 = 1e4)
%     n1 = 100;
%     g1 = linspace(0,1,n1);
%     % redistribute along string
%     dx2 = diff(x,1,2);
% 
%     dd = sum(dx2.^2);
%     dd = sqrt([0 dd]);
% 
%     ll = cumsum(dd);
%     ll = ll/ll(end);
%     x = interp1q(ll',x',g1')';
    
    dx = (x(:,3:n1) - x(:,1:(n1-2)))/2;
    dW = zeros(p, n1);
    conv_flag = zeros(n1,1);
    
    r = zeros(n, n1);
    
    H_arr = zeros(1, n1);
    for j = 2:(n1-1)
    %parfor (j = 2:(n1-1), 3)
        dg = deterministic_force(system, x(:, j));
        dp = dot(dg, dx(:, j-1));
        % figure out which direction we should compute the quasipotential
        % in
        if ~isempty(guess)
            subguess = guess(:, j-1);
        else
            subguess = -dg;
        end
        [dW(:, j), conv_flag(j)] = dW_fun(system, x(:,j), -sign(dp)*dx(:,j-1), subguess);
        [H_arr(j), r(:, j)] = H(dW(:,j), system.S(x(:,j)), system.V(x(:,j)));
    end
    % fill in dx at endpoints
    dx = [x(:,2)-x(:,1) dx x(:,end)-x(:,(end-1))];
    W = cumsum(sum(0.5*(dW(:, 2:end).*dx(:, 2:end) + dW(:, 1:(end-1)).*dx(:,1:(end-1))),1));
    W = [0 W];
end