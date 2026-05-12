% Produces the switching matrix given bead positions
% Right now just hardcode the 3-bead case
function [S, r] = switching(x)
    x = reshape(x, 3, []);
    [d, p] = size(x);
    %a = @(d) 2*(abs(d)<0.75);
    %a = @(d) 2*max(0.75-abs(d), 0);
    %a = @(d) 0;
    %a = @(d) 2*exp(-(abs(d)).^2);
    a = @(d) 2./(1+exp(20*(abs(d)-0.75)));
    c = 0.5;
    
    configs = state_configurations(p);
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