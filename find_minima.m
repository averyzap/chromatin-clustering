function out = find_minima(system, start, verbose, fix)
    d = system.d;
    x = start;
    dt = 2e-1;
    maxsteps = 10000;
    
    if nargin < 3
        verbose = 0;
    end
    
    tol = 1e-12;
    
    %start = start + 2*tol*normrnd(0,1,size(start));
    
    if nargin < 4
        fix  = 1;
    end
    
    for i=1:maxsteps
        xold = x;
        x = x + dt*deterministic_force(system, x);
        
        if fix
            if d==2
                x(5) = 0;
            end
        end
        
        if max(abs(xold-x)) < tol
            break
        end
        
        if verbose && mod(i,100)==0
            fprintf("Iteration: %d   change: %e\n",i,max(max(abs(x-xold))));
        end
    end
    
    out = x;
end