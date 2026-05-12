function [out] = M(w, S, V, shift)
    if nargin < 4
        shift = 0;
    end

    w = reshape(w, 1, []);
    
    n = length(S);
    shift_mat = -n*shift*eye(n)+shift;
    %shift_mat = shift*(tril(S,-1)~=0);
    
    out = S + diag(w*V) + (w*w')*eye(n) + shift_mat;
end