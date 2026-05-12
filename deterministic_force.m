function out = deterministic_force(system, x)
%     if nargin < 2
%         d = 3;
%     end
    V = drift_mat(system, x);
    S = switching(system, x);
    if any(any(isnan(S) | isinf(S)))
        x=1;
    end
    r = null(S);
    r = r/sum(r);
    out = reshape((V*r), size(x));
end