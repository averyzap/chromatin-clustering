function out = find_minima(system, start)
    x = start;
    dt = 1e-1;
    maxsteps = 10000;
    
    tol = 1e-12;
    
    %start = start + 2*tol*normrnd(0,1,size(start));
    
    for i=1:maxsteps
        xold = x;
        x = x + dt*deterministic_force(system, x);
        
        if max(abs(xold-x)) < tol
            break
        end
        
        if mod(i,100)==0
            fprintf("Iteration: %d   change: %e\n",i,max(max(abs(x-xold))));
        end
    end
    
    out = x;
end