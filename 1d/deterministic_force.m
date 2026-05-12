function out = deterministic_force(system, x)
    V = system.V(x);
    [~,r] = system.S(x);

    %r = null(S);
    %r = r/sum(r);
    out = reshape((V*r), size(x));
end