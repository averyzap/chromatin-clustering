% Produces the switching matrix given bead positions
function [S, r] = switching(system, x, configs)
    d = system.d;
    a = system.a;
    
    x = reshape(x, d, []);
    [d, p] = size(x);
    c = 0.5;
    
    if nargin < 3
        sc = memoize(@state_configurations);
        configs = sc(length(x));
    end
    n = size(configs,2);
    S = zeros(n);
    
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
                if s2(f(1))==f(1)
                    % bond broke from s1 to s2: this happens at rate given by c
                    S(i, j) = c;
                else
                    % bond formed - rate depends on distance
                    S(i, j) = a(norm(x(:, f(1))-x(:, f(2))));
                end
            end
        end
    end
    
    S = S - diag(sum(S,2));
    S = S';
    
    if nargout == 2
        r = null(S);
        r = r/sum(r);
    end
end