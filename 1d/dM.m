% derivative of M wrt i'th momentum variable
function [dMv, dMv2] = dM(i, w, S, V)
    dMv = diag(V(i, :)) + 2*w(i)*eye(size(S,1));
    if nargout == 2
        dMv2 = 2*eye(size(S,1));
    end
end