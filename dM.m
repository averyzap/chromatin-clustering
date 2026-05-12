% derivative of M wrt all momentum variable
function [dMv, dMv2] = dM(system, x, w, S, V, with_x)
    p = length(w);
    
    if nargin < 6
        with_x = 0;
        imax = p;
    else
        imax = 2*p;
    end

    n = size(S, 1);
    dMv = zeros(n,n,imax);
    if nargout == 2
        dMv2 = zeros(size(dMv));
    end
    
    for i=1:p
        dMv(:,:,i) = diag(V(i, :)) + 2*w(i)*eye(size(S,1));
        if nargout == 2
            dMv2(:, :, i) = 2*eye(size(S,1));
        end
    end
       
    if with_x
        % need derivatives from S and A
        [dSm] = dS(system, x);
        dMv(:,:,(p+1):end) = dSm;
        dVm = zeros([size(V) p]); 

        % TODO: do A derivative numerically
        for xi = 1:p
            xv = x;
            xv(xi) = xv(xi) + 1e-5;
            Vp1 = drift_mat(xv, system.d);
            Sp1 = switching(system, xv);
            xv(xi) = xv(xi) - 2e-5;
            Vp2 = drift_mat(xv, system.d);
            Sp2 = switching(system, xv);
            tmp = (Vp1-Vp2)/(2e-5);
            tmp2 = (Vp1 - 2*V + Vp2)/(1e-10);
            asdf=1;
           %  dVm(:,:,xi)
            dMv(:,:,p+xi) = dMv(:,:,p+xi) + diag(w' * tmp);
            dMv2(:,:,p+xi) = dMv2(:,:,p+xi) + diag(w' * tmp2);

            Stmp2 = (Sp1 - 2*S + Sp2)/(1e-10);
            dMv2(:,:,p+xi) =  dMv2(:,:,p+xi) + Stmp2;
        end
    end
end