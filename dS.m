function [out] = dS(system, x)
    d = system.d;
    
    if isfield(system, 'a')
        if ~isfield(system, 'da')
            error('Must provide derivative of affinity function in system');
        end
        da = system.da;
    else
        da = @(x) -40*x*exp(5*(4*abs(x)+3))/(abs(x)*(exp(20*abs(x))+exp(15))^2);
    end
    
    x = reshape(x, d, []);
    [d, p] = size(x);
    %a = @(d) 2*(abs(d)<0.75);
    %a = @(d) 2*max(0.75-abs(d), 0);
    %a = @(d) 0;
    %a = @(d) 2*exp(-3*(abs(d)).^2);
    %a = @(x) 2./(1+exp(20*(abs(x)-0.75)));
    %a = @(d) 4./(1+exp(20*(abs(d)-0.75)));
    c = 0.5;
    
    configs = state_configurations(p);
    n = size(configs,2);
    S = zeros(n);
    out = zeros(n,n,d, p);
    
    %v2 = rem(xi-1, 2)+1;
    %v1 = (xi-v2)/2 + 1;
    
    for i=1:n
        for j=1:n
            if i==j
                continue
            end
            s1 = configs(:, i);
            s2 = configs(:, j);
            
            changed = (s1~=s2);
            f = find(changed);
            
            if length(f)==2
                if ~(s2(f(1))==f(1))
                    
                    dv = norm(x(:, f(1))-x(:, f(2)));
                    dav = da(dv);
                    
                    % loop through all relevant coordinates and add them to
                    % the derivative
                    for kk = 1:d
                        du = 1/(2*dv) * 2 *(x(kk, f(1))-x(kk, f(2)));
                        
                        out(i, j, kk, f(1)) = out(i, j, kk, f(1)) + dav*du;
                        out(i, j, kk, f(2)) = out(i, j, kk, f(2)) - dav*du;
                    end
                end
            end
        end
    end
    
    out = reshape(out, n,n,[]);
    for i=1:(d*p)
        out(:,:,i) = (out(:,:,i) - diag(sum(out(:,:,i), 2)))';
        %S = S - diag(sum(S,2));
    end
    %S = S - diag(sum(S,2));
    %out = S';
end