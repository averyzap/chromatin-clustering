function [out, r] = H(w, S, V, n, shift, filter)
    if nargin < 4
        n = 1;
    end
    if nargin < 5
        shift = 0;
    end
    
    Mv = M(w, S, V, shift);
    
    if nargin < 6
        filter = ones(length(Mv), 1);
    end
    
    filter = logical(filter);
    
    Mv = Mv(filter, :);
    Mv = Mv(:, filter);
    
    [v,d] = eig(Mv);
    e = diag(d);
    %
    
    % maxk missing workaround
    if n == 1
        [out, I] = max(real(e));
    else
        out = maxk(e, n,'ComparisonMethod','real');
    end
    
    if nargout == 2
        r = v(:, I);
        r = r/sum(r);
    end
    
    if ~isreal(out)
        warning('Imaginary eigenvalues');
        out = real(out);
    end
end