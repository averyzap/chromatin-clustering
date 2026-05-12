function V = well_drift_mat(x)
    % TODO: extend to higher dimensionality
    assert(isscalar(x));
    % two switching states - well on vs off
    V = zeros([1 2]);
    
    % binding spring force
    kv = 5;
    % excluded volume strength
    a_ev = 3;
    % excluded volume distance
    c_ev = 0.5;
    
    phi = 1;
    
    
    V(1) =  a_ev * (x)...
            * exp(-norm(x).^2 / c_ev);
    V(2) = V(1) - kv*x;
    
    V = V/phi;
end