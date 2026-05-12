function [out, hess] = gradH(system, x, w, S, V, k, with_x)
    if nargin < 6
        k = 1;
    end
    if nargin < 7
        with_x = 0;
        % only use 3 arguments
        V = w;
        S = x;
        w = system;
    end
    
    
    Mvb = M(w,S,V);
    
    %if nargin < 6
        filter = ones(length(Mvb), 1);
    %end
    
    filter = logical(filter);
    
    Mvb = Mvb(filter, :);
    Mvb = Mvb(:, filter);
    
    n = length(S);
    p = length(w);
    
    if with_x
        imax = 2*p;
    else
        imax = p;
    end
    
    out = zeros(imax, k);
    hess = zeros(imax, imax, k);
    
    
    [v, d] = eig(Mvb);
    
    [vt, dt] = eig(Mvb');
    if k==1
        [lv, I] = max(real(diag(d)));
        [lvt, It] = max(real(diag(dt)));
    else
        [lv, I] = maxk(real(diag(d)),k);
        [lvt, It] = maxk(real(diag(dt)), k);
    end
    
    [dMv_arr, dMv2_arr] = dM(system, x, w, S, V, with_x);

    for lvi = 1:k
        l = lv(lvi);
        Mv = Mvb-l*eye(size(Mvb));
        Mvi = pinv(-Mv);
        u0 = v(:, I(lvi));
        v0 = vt(:, It(lvi));
        
        if dot(u0, v0) == 0
            % failure - just fall back to numerical for now
            warning('Analytic gradient failure');
            [out, hess] = test_hess(w, S, V, k, shift);
            return;
        end

        if nargout == 2
            %hess = zeros(size(Mv));
            try
                K0 = eye(size(Mv)) - (u0*v0')/(v0'*u0);
            catch
                error('asdf');
                x=1;
            end
            %dMv2 = 2*eye(size(Mv));
            
        end

        for i=1:imax
            [dMv] = dMv_arr(:,:,i);
            
            dMv = dMv(filter, :);
            dMv = dMv(:, filter);
            
            out(i, lvi) = (v0' * dMv *u0)/(v0'*u0);
        end

        if nargout == 2
            for i=1:imax
                for j=i:imax
                    hess(i, j, lvi) = 2*(v0' * dMv_arr(:, :, i) * K0 * Mvi * (K0 * dMv_arr(:, :, j) * u0))/(v0'*u0);
                    if i~=j
                        hess(j, i, lvi) = hess(i, j, lvi);
                    end
                end
                % extra term on diagonal
                hess(i, i, lvi) = hess(i, i, lvi) + (v0' * dMv2_arr(:,:,i) *u0)/(v0'*u0);
            end
        end
    end
    
    if any(any(isnan(out))) || any(any(any(isnan(hess)))) 
        warning('Analytic gradient failure');
        [out, hess] = test_hess(w, S, V, k, shift);
    end
    
    if ~isreal(out) || ~isreal(hess)
        out = real(out);
        hess = real(hess);
    end
end