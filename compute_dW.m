function [dW, flag, iter, Hppv] = compute_dW(system, x, dx, guess, verbose)
    ndx = norm(dx);
    udx = dx / ndx;
    d = system.d;
    sc = memoize(@state_configurations);
    configs = sc(length(x)/d);
    [S,r] = switching(system, x, configs);
    V = drift_mat(system, x, configs);
    % iterative update procedure as per Jay's paper
    maxiters = 1000;
    tol = 1e-12;
    H_tol = 1e-6;
    angle_tol = 1e-5;
    gamma = 1;
    
    if nargin < 5
        verbose = 0;
    end
    
    flag = 1;
    
    dt = 9000;
    
    change = nan;
    change2 = nan;
    
    delta_arr = zeros(1, maxiters);
    
    Hguess = H(guess, S, V);
    [guess, Hguess, Hpguess] = rescue(guess, Hguess, H_tol, S, V);
    
    p = guess;
    pold = guess;
    
    % backup guess from dx
    gv = 5;
    Hv = -1;
    while Hv < 1
        pv = gv * dx;
        Hv = H(pv, S, V);
        gv = gv*2;
    end
    [pv, Hv, Hp] = rescue(pv, Hv, H_tol, S, V);
    % figure out which guess is better
    og = dot(guess, dx);
    ov = dot(pv, dx);
    if ov > og
        obj = ov;
        p = pv;
        Hpnew = Hp;
    else
        obj = og;
        p = guess;
        Hpnew = Hpguess;
    end
    angle = acos(dot(udx, Hpnew)/norm(Hpnew));
    
    for iter = 1:maxiters
        % try a Newton step from our current best point
        Hv = H(p, S, V, 1);
        
        % Compute gradient and hessian of hamiltonian
        [Hp, Hpp] = gradH(p, S, V, 1);
        
        % compute Lagrange multiplier
        a = (dot(Hp, Hpp \ Hp) - 2*Hv)/dot(dx, Hpp \ dx);
        if a < 0
            a = 0;
        end
        lambda = ndx * sqrt(a);
        
        
        % attempt update per Newton method
        pnew = p + gamma*(Hpp \ (lambda * udx - Hp));
        
        Hnew =  H(pnew, S, V);
        
        [pnew, Hnew, Hpnew] = rescue(pnew, Hnew, H_tol, S, V);

        
        change = max(abs(pnew-p));
        
        if change < tol*gamma || (angle < angle_tol && abs(Hnew) < H_tol)
            dW = pnew;
            flag = 0;
            
            if angle > angle_tol
                warning('high angle at convergence')
                dW = zeros(size(dW));
                flag = 1;
            end
            if verbose
                fprintf("Iter: %3d\tH:%e\tdx*p:%.16e\tAngle:%.3e\tdelta:%e\td2:%e\tgamma:%e\tl:%.2e\n",iter, Hnew, dot(dx, pnew), angle, change, change2, gamma, lambda);
            end
            if nargout >= 4
                [~, Hpp] = gradH(dW, S, V);
                %[~, Hpp] = test_hess(dW, S, V, 1, 1e-5);
                Hppv = norm(Hpp);
            end
            return
        end
        
        % check if we need to do a fallback
        cur_obj = dot(dx, pnew);
        
        % don't let dt get stuck too small
        if dt < 1e-8
            dt = 1e-6;
        elseif dt < 1e+4
            dt = 1.4*dt;
        end
        %dt = 5;
        siter = 0;
        while cur_obj < obj
            % it's not an improvement - fallback
            dt = 0.75*dt;
            % perhaps adjusting direction would help
            v = dx - dot(dx, Hp)/(norm(Hp)^2)*Hp;
            v = v/norm(v);
            if change2 < change
                % we're going back and forth, try to correct
                v2 = p - pold;
                v2 = v2/norm(v2);
                v = 0.5*(v+v2);
                v = v/norm(v);
            end
            pnew = p + dt*v;
            Hnew = H(pnew, S, V);
            [pnew, Hnew, Hpnew] = rescue(pnew, Hnew, H_tol, S, V);
            % fallback should be an improvement
            cur_obj = dot(dx, pnew);
            if siter > 100
                asdf = pnew;
                dW = 0;
                flag = 1;
                return;
            end
            siter = siter+1;
        end
        
        if verbose && (siter > 5)
            disp(siter);
        end
        
        change = max(abs(pnew-p));
        change2 = max(abs(pnew-pold));
        
        delta_arr(iter) = change;
        angle = acos(dot(udx, Hpnew)/norm(Hpnew));
        
        if verbose
            fprintf("Iter: %3d\tH:%.3e\tdt:%.3e\tdx*p:%.8e\tAngle:%.3e\tdelta:%e\td2:%e\tgamma:%e\tl:%.2e\n",iter, Hnew, dt, dot(dx, pnew), angle, change, change2, gamma, lambda);
        end
        pold = p;
        p = pnew;
        obj = cur_obj;
        
        %disp(p')
        
    end
    
    if iter == maxiters
        if angle < angle_tol
            if verbose
                warning('Max iterations but low angle');
            end
            flag=0;
        else
            p = zeros(size(p));
            if verbose
                warning('Iteration did not converge')
            end
            flag=1;
        end
    end
    
    dW = p;
    [~, Hpp] = gradH(p, S, V);
    Hppv = norm(Hpp);
end

% Rescue operation
function [pnew, Hnew, Hpnew] = rescue(pnew, Hnew, H_tol, S, V)
    Hpnew = [];
    iter = 0;
    while abs(Hnew) > H_tol || Hnew > 1e-13
        if iter>100
            asdf=1;
        end
        
        Hpnew = gradH(pnew, S, V);
        delta = Hpnew*Hnew/dot(Hpnew,Hpnew);
        if Hnew > 0 && abs(Hnew) < H_tol
            delta = 1.1*delta;
        end
        pnew = pnew - delta;
        Hnew = H(pnew, S, V);
%         if Hnew > 0 && abs(Hnew) < H_tol
%             pnew = pnew - H_tol*Hpnew;
%             Hnew = H(pnew, S, V);
%         end
        iter = iter+1;
    end
    
    % we're expected to still return this even if we didn't enter the
    % above loop
    if isempty(Hpnew)
        Hpnew = gradH(pnew, S, V);
    end
end