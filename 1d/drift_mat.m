% Returns the matrix V that describes state-dependent drift
function V = drift_mat(x, d)
    if nargin < 2
        d = 3;
    end
    x = reshape(x, d, []);
    % binding spring force
    kv = 5;
    % confinement spring force
    kc = 0.1;
    % excluded volume strength
    a_ev = 3;
    % excluded volume distance
    c_ev = 0.5;
    phi=0.8;
    configs = state_configurations(length(x));
    n = size(configs,2);
    p = size(configs,1);
    V = zeros(d, p, n);
    for i=1:n
        state = configs(:, i);
        for k=1:p
            V(:, k, i) = kv*(x(:, state(k))-x(:, k)) - kc*x(:,k)*norm(x(:,k))^2;
            % add in excluded volume
            for l=1:p
                if l==k
                    continue
                end
                V(:,k, i) = V(:,k, i) + a_ev * (x(:,k) - x(:,l))...
                    * exp(-norm(x(:,k)-x(:, l)).^2 / c_ev);
            end
            V(:,k, i) = V(:,k, i)/phi;
        end
    end
    V = reshape(V, [d*p n]);
end