function outstruct = string_v2(x, dt, n1, dW)
if ndims(x) == 3
    x = reshape(x, [], size(x, 3));
end
Np = size(x,1);
% number of images along the string (try from  n1 = 3 up to n1 = 1e4)
%n1 = 30;
g1 = linspace(0,1,n1);
% redistribute along string
dx = diff(x,1,2);

dd = sum(dx.^2);
dd = sqrt([0 dd]);

ll = cumsum(dd);
ll = ll/ll(end);
x = interp1q(ll',x',g1')';

xOld = x;

x_arr = {};

if any(isnan(x))
    x=1;
end


% number of steps of steepest descent
nstepmax = 4000;

plot_during=0;
make_movie = 1;
nn = 5;   % when to plot

if make_movie
    v = VideoWriter('string_evolution.mp4', 'MPEG-4');
    open(v);
%     if ~plot_during
%         f = figure('visible','off');
%     end
    frames = [];
    bounds = [-1 1 -1 1 -1 1];
end

% time-step (limited by the ODE step but independent of n1)
%dt = 0.01;
dt_delta = 0.1;

% parameter used as stopping criterion (not always used)
diff_tol1 = 1e-4;
exitflag = 0;   % changes to 1 if convergence reached

% Initialize dW to be dx (used as initial guess in recomputation of dW
% later
%dW = (x(:,3:n1) - x(:,1:(n1-2)))/2;
%dW2 = -(x(:,3:n1) - x(:,1:(n1-2)))/2;


if nargin < 4
    dW = zeros(Np, n1);
    
%     for j=2:(n1-1)
%         dW(:, j) = deterministic_gradient(x(:, j))/10;
%     end
else
    dW = interp1q(ll',dW',g1')';
    %dW = dW(:, 2:(end-1));
end
% approximation of the Hessian
d2W = zeros(size(dW));
    
diffX = zeros(size(x));

W_arr = zeros(1, nstepmax);
err_arr = zeros(1, nstepmax);

dt_gmam = dt;

errX = 1;

gmam = 0;

climbing = 1;

fix_beads = 1;

if climbing
    end_val = n1;
else
    end_val = n1-1;
end

gmam_start = 20;

if ~climbing
    dW(:, end) = 0;
end

for nstep = 1:nstepmax
    x0 = x;
    
    if nstep <= gmam_start
        dW_fun = @compute_dW_v4;
    else
        dW_fun = @compute_dW_v4;
    end
    
    dphi = (x(:,2:n1) - x(:,1:end-1))/2;
    dphi = [dphi x(:, end)-x(:,end-1)];
    %dphi = x(:, 2:n1) - x(:, 1:(n1-1));
    
    dx = (x(:,3:n1) - x(:,1:(n1-2)))/2;
    % fill in dx at endpoints
    dx = [x(:,2)-x(:,1) dx x(:,end)-x(:,(end-1))];
    
    % gradient descent for each image except endpoints
    conv_flag = zeros(n1-2, 1);
    dWOld = dW;
    dp_arr = zeros(n1, 1);
    %parfor (j = 2:(n1-1), 2)
    for j=2:end_val
        if nstep == 1 && j > 2
            guess = dW(:, j-1);
        else
            guess = dW(:, j);
        end
        
        % figure out which direction we should compute the quasipotential
        % in
        dg = deterministic_gradient(x(:, j));
        dp = dot(dg, dphi(:, j-1));
        dp_arr(j) = dp;
        %dW2(:, j) = compute_dW(x(:,j), -sign(dp)*dx(:,j-1), dg/2);
    
        %[dW(:, j-1), conv_flag(j-1)] = compute_dW(x(:,j),  dphi(:,j-1),  guess);
        [tmp, conv_flag(j-1)]=  compute_dW_v6(x(:,j),  -sign(dp)*dphi(:,j-1),  guess);
        %[tmp, conv_flag(j-1)] = compute_dW_v4(x(:,j),  dphi(:,j-1),  guess);
        if conv_flag(j-1)==0
            dW(:, j) = tmp;
        else
            if nstep < 5
                [tmp2, conv2] = compute_dW(x(:,j),  -sign(dp)*dphi(:,j-1),  guess);
            else
            conv2=1;
            end
            if conv2 == 0
                dW(:, j) = tmp2;
            else
                fprintf("Failed iteration %d\n",j);
                dW(:, j) = 0;
            end
        end
%         elseif conv_flag(j-1)==3
%             % blew up, generally because eigenvalues not differentiable
%             fprintf("Failed iteration %d\n",j);
%             dW(:, j) = compute_dW(x(:,j),  -sign(dp)*dphi(:,j-1),  guess);
%         else
%             fprintf("Failed iteration %d\n",j);
%             dW(:, j) = compute_dW(x(:,j),  -sign(dp)*dphi(:,j-1),  guess);
%             %dW(:, j-1) = 0.5*(dWOld(:, j-2) + dWOld(:, j));
% %             fprintf("Didn't converge on image %d\n",j);
% %             [dW(:, j-1), conv_flag(j-1)] = compute_dW(x(:,j),  -sign(dp)*dphi(:,j-1),  guess);
% %             if conv_flag(j-1)
% %                 a=1;
% %                 dW(:, j-1) = 0.5*(dWOld(:, j-1) + dWOld(:, j+1));
% %             end
%         end
        
        % compute the estimate of the hessian
        eps=5e-3;
        if gmam && nstep >= gmam_start
            [d2W(:, j), d2w_flag] = compute_dW_v6(x(:, j) + eps*dW(:, j), -sign(dp)*dphi(:,j), dW(:, j));
            d2W(:, j) = (d2W(:, j) - dW(:, j))/eps;
            
            % sometimes it seems like this doesn't get computed well
%             d2W_diff = norm(d2W(:, j-1), Inf)/norm(dW(:, j-1), Inf);
%             if 0 && (d2W_diff > 5 || d2W_diff == 0)
%                 d2W(:, j-1) = dW(:, j-1);
%             else
%                 a=1;
%             end
            if flag ~= 0
                a=1;
            end
        end
        
        %dW2(:, j-1) = compute_dW(x(:,j), -dphi(:,j-1), dW2(:,j-1));
    end
    
    % test debug - fit a convex shell over the values. first find where
%     vn_arr = [0 vecnorm(dW) 0];
%     switch_point = max(find(dp_arr<0));
%     ind = (1:switch_point)';
%     P = [ind vn_arr(ind)'];
%     k = convhull(P);
%     k = flip(k(find(k==switch_point,1):end));
%     
%     dWs = [zeros(9,1) dW];
%     dWs = dWs(:, k);
%     xs = x(:, k);
%     g2 = linspace(0,1,switch_point);
%     % redistribute along string
%     dxs = diff(xs,1,2);
% 
%     dd = sum(dxs.^2);
%     dd = sqrt([0 dd]);
% 
%     ll = cumsum(dd);
%     ll = ll/ll(end);
%     dWs = interp1q(ll',dWs',g2')';
%     
%     dW(:,1:(switch_point-1)) = dWs(:, 2:end);
    
    %W = cumsum(sum(dW.*dx));
%     dWp = [zeros(9,1) dW zeros(9,1)];
%     W = cumsum(sum(0.5*(dWp(:, 2:end).*dx(:, 2:end) + dWp(:, 1:(end-1)).*dx(:,1:(end-1)))));
%     W = [0 W];
%     Wmax = max(W);

    
    if gmam && nstep >= gmam_start
        dt_gmam = dt;
        
        % GMAM formula update
        dx = x(:,2:n1) - x(:,1:n1-1);
        Lk = sum(vecnorm(dx));

        gk = vecnorm(dW)/Lk;
        %gk = [0 gk 0];
        
        % Compute the RHS vector
        d2Wp = d2W(:, 2:end-1);
        b = -d2Wp + gk(2:(end-1)).*(gk(3:end)-gk(2:(end-1))).*dx(:, 2:end)./(2/n1^2) ...
            + gk(2:(end-1)).*(gk(2:(end-1))-gk(1:(end-2))).*dx(:, 1:(end-1))./(2/n1^2) ...
            + x(:, 2:(n1-1))/dt_gmam;
        
        gk2 = gk(2:(end-1)).^2;
        % Construct the LHS matrix
        A = (1/n1^2) * (diag(gk2(1:(end-1)).*ones(1,n1-3), 1) + diag(gk2(2:end).*ones(1,n1-3),-1) - 2*diag(gk2.*ones(1,n1-2)));
        % duplicate A over all 9 dimensions
        tmp = repmat({A},1,Np);
        A = blkdiag(tmp{:});
        A = -A + (1/dt_gmam)*eye(size(A));
        
        % reshape A and b to facilitate matrix problem
        %A = reshape(A, (n1-2)*Np, (n1-2)*Np);
        b = reshape(b',[],1);

        xnew = reshape(A \ b, [], 9)';
        xnew2 =  x(:,2:(n1-1)) - dt*(dW(2:n1-1));
        x(:, 2:(end-1)) = xnew;
        
        % for climbing, handle the final image like the string method
        if climbing
            % change final image to have the climbing direction
            last_change = -dt*dW(:,end);
            last_change = last_change - 1.5*dphi(:,end)*dot(last_change, dphi(:, end))/(norm(dphi(:,end))^2);
            
            x(:, end) = x(:, end) + last_change;
        end
    else
        % string method update
        %dt_sub = min(dt, sqrt(errX));
        change = -dt*dW;
        
        % debug: try removing parallel component from everything
%         for j=2:(n1-1)
%             change(:, j) = change(:, j) - dphi(:,j-1)*dot(change(:,j), dphi(:, j-1))/(norm(dphi(:,j-1))^2);
%         end
        
        if climbing
            % change final image to have the climbing direction
            nold = norm(change(:,end));
            change(:, end) = change(:, end) - 1.5*dphi(:,end)*dot(change(:,end), dphi(:, end))/(norm(dphi(:,end))^2);
            assert(norm(change(:,end))<=nold);
        else
            % make sure we're not screwing this up
            assert(all(change(:, 1)==0) && all(change(:,end)==0));
        end
        
        if fix_beads
            change(4:end,:)=0;
        end
        
        x = x + change;
    end
    
    Wmax = getW(x,dW);
        
    
    % redistribute along string
    dx = x(:,2:n1) - x(:,1:n1-1);
    
    dd = sum(dx.^2);
    dd = sqrt([0 dd]);
    
    ll = cumsum(dd);
    string_len = ll(end);
    ll = ll/ll(n1);
    xnew = interp1q(ll',x',g1')';
    x=xnew;
    
    
    % --- exit condition
    diffold = diffX;
    diffX = x-x0;
    diffX_arr = diffX(:);
    diffFromStart = x-xOld;
%     
%     if nstep > 95
%         plot(vecnorm(diffX));
%         pause;
%     end

    
    
    errX = norm(diffX(:), Inf);
    errX2 = norm(diffX(:), 1);
    errFromStart = norm(diffFromStart(:), Inf);
    errEnd = norm(diffX(:, end), Inf);
    err_align = dot(diffold(:), diffX(:))/(norm(diffold(:))*norm(diffX(:)));
    fprintf("Iteration: %d\tW max:%e\tError:%e\tAlignment:%f\tLength: %e\n", nstep,Wmax, errX, err_align, string_len);
    W_arr(nstep) = Wmax;
    err_arr(nstep) = errX;
    
    x_arr{numel(x_arr)+1} = x;
    
    if nstep > 10 && err_align < 0.99
        dt = 1e-3;
    end
    
    if nstep > 500 && err_align < 0
        a=1;
    end
    
%     if nstep > 100 && err_arr(nstep) > err_arr(nstep-5)
%         W_arr = W_arr(1:nstep);
%         err_arr = err_arr(1:nstep);
%         break;
%     end

    if nstep == 500
        xOld = x;
    end
    
    if (make_movie || plot_during) && mod(nstep, nn) == 0
        plot3dbeads(x);
        axis(bounds);
        %test_eigvals(x,dW);
        F = getframe;
        writeVideo(v, F);
        
        if plot_during
            % pause to allow viewing of the current state
            pause;
        end
    end
    
    %errdV = max(abs(dW));
    if (nstep > gmam_start && err_align < 0)
       %dt = dt/2;
    elseif nstep > gmam_start && err_align > 0.99
        %dt = 1.1*dt;
    end
    if errX < diff_tol1*dt 
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

dx = x(:,2:n1) - x(:,1:n1-1);

dd = sum(dx.^2);
dd = sqrt([0 dd]);

%x = reshape(x, 3, [], n1);

%[x,dW,exitflag, W_arr, err_arr, x_arr]
outstruct = [];
outstruct.x = x;
outstruct.dW = dW;
outstruct.exitflag = exitflag;
outstruct.W_arr = W_arr;
outstruct.err_arr = err_arr;
outstruct.x_arr = x_arr;

end

function Wmax = getW(x, dW)
    dx = (x(:,3:end) - x(:,1:(end-2)))/2;
    % fill in dx at endpoints
    dx = [x(:,2)-x(:,1) dx x(:,end)-x(:,(end-1))];
    %W = cumsum(sum(dW.*dx));
    dWp = dW;
    W = cumsum(sum(0.5*(dWp(:, 2:end).*dx(:, 2:end) + dWp(:, 1:(end-1)).*dx(:,1:(end-1)))));
    W = [0 W];
     Wmax = max(W);
end