function [S, r] = well_switching(x, a)
    if nargin < 2
        a = @(d) 2*exp(-3*d.^2);
        %a = @(d) 2./(1+exp(20*(abs(d)-0.75)));
        %a = @(d) 4./(1+exp(20*(abs(d)-0.75)));
        %a = @(d) 2;
    end
    c = 0.5;
    
    % binding
    S(1,2) = a(x);
    S(2,1) = c;
    
    S = S - diag(sum(S,2));
    S = S';
    
    if nargout == 2
        r = null(S);
        r = r/sum(r);
    end
end