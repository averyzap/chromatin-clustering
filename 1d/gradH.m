function [out, hess] = gradH(w, S, V, k, shift)
    if nargin < 4
        k = 1;
    end
    if nargin < 5
        shift = 0;
    end
    Mvb = M(w,S,V, shift);
    n = length(S);
    p = length(w);
    out = zeros(p, k);
    hess = zeros(p, p, k);
    
    
    [v, d] = eig(Mvb);
    
    [vt, dt] = eig(Mvb');
    if k==1
        [lv, I] = max(real(diag(d)));
        [lvt, It] = max(real(diag(dt)));
    else
        [lv, I] = maxk(real(diag(d)),k);
        [lvt, It] = maxk(real(diag(dt)), k);
    end
    
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
                x=1;
            end
            dMv2 = 2*eye(size(Mv));
            dMv_arr = zeros([size(Mv) p]);
        end

        for i=1:p
            [dMv] = dM(i,w,S,V);
            out(i, lvi) = (v0' * dMv *u0)/(v0'*u0);

            % store differentials for later use in Hessian
            if nargout == 2
                dMv_arr(:, :, i) = dMv;
            end
        end

        if nargout == 2
            for i=1:p
                for j=i:p
                    hess(i, j, lvi) = 2*(v0' * dMv_arr(:, :, i) * K0 * Mvi * (K0 * dMv_arr(:, :, j) * u0))/(v0'*u0);
                    if i~=j
                        hess(j, i, lvi) = hess(i, j, lvi);
                    end
                end
                % extra term on diagonal
                hess(i, i, lvi) = hess(i, i) + (v0' * dMv2 *u0)/(v0'*u0);
            end
        end
    end
    
    if any(any(isnan(out))) || any(any(any(isnan(hess)))) 
        [out, hess] = test_hess(w, S, V, k, shift);
    end
    
    if ~isreal(out) || ~isreal(hess)
        out = real(out);
        hess = real(hess);
    end
end