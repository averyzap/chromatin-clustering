% Returns the matrix V that describes state-dependent drift
function [V, dVdx] = drift_mat(system, x, configs)
    d = system.d;
    
    x = reshape(x, d, []);
    [d, p] = size(x);
    % binding spring force
    kv = system.kv;
    % confinement spring force
    kc = system.kc;
    % excluded volume strength
    a_ev = system.a_ev;
    % excluded volume distance
    c_ev = system.c_ev;
    phi = system.phi;
    if nargin < 3
        sc = memoize(@state_configurations);
        configs = sc(length(x));
    end
    n = size(configs,2);
    p = size(configs,1);
    V = zeros(d, p, n);
    if nargout == 2
        dVdx = zeros([size(V) numel(x)]);
    end
    for i=1:n
        state = configs(:, i);
        for k=1:p
            V(:, k, i) = kv*(x(:, state(k))-x(:, k)) - kc*x(:,k)*norm(x(:,k))^2;
            if nargout == 2
                if state(k) ~= k
                    for jj = 1:d
                        v1 = (k-1)*d + jj;
                        dVdx(jj, k, i, v1) = dVdx(jj, k, i, v1) - kv;
                    end
                end
            end
            % add in excluded volume
            for l=1:p
                if l==k
                    continue
                end
                V(:,k, i) = V(:,k, i) + a_ev * (x(:,k) - x(:,l))...
                    * exp(-norm(x(:,k)-x(:, l)).^2 / c_ev);
                if nargout == 2
                    for jj = 1:d
                        v1 = (k-1)*d + jj;
                        dVdx(jj, k, i, v1) = dVdx(jj, k, i, v1) + ...
                            a_ev *  exp(-norm(x(:,k)-x(:, l)).^2 / c_ev) + ...
                            a_ev *  (x(:,k) - x(:,l)) * exp(-norm(x(:,k)-x(:, l)).^2 / c_ev) * 1 ... 
                                % TODO: finish this
                    end
                end
            end
            V(:,k, i) = V(:,k, i)/phi;
        end
    end
    V = reshape(V, [d*p n]);
end