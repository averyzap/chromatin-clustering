% Right now assumes everything is one-dimensional
% TODO: let it handle transitions within less than a timestep
function [pos_t, state_t, step] = simulate_system(system, x, dt, nt, eps, stop)
    %rng(356);
    d = system.d;
    x = reshape(x, d, []);
    p = size(x,2);
    nb = p*(p-1)/2;
    
    if nargin < 6
        stop = [];
    end
    
    buffer_size = 1000000;
    
    fix_beads = 0;
    
    % binding spring force
    kv = system.kv;
    % confinement spring force
    kc = system.kc;
    % excluded volume strength
    a_ev = system.a_ev;
    % excluded volume distance
    c_ev = system.c_ev;
    phi = system.phi;
    
    if isfield(system, 'a')
        a = system.a;
    else
        error('System is missing affinity function.')
        %a = @(x) 2./(1+exp(20*(abs(x)-0.75)));
    end
    c = 0.5;
    
    t = 0;
    
    old_bonds = 1;
    
    break_mean = eps/c;
    
    planestop = isfield(stop, 'planes');
    
    if fix_beads
        %dold = norm(x(1:2)-x(3:4));
        %dold = norm(x(:,1);
        % rotate system to maintain symmetry fixing
        xv = x(5:6);
        %target = [0 -1];
        %theta = -acos(dot(xv,target)/norm(xv));
        r = norm(xv);
        c = xv(2)/r;
        s = xv(1)/r;
        R = [c -s; s c];
        x(1:2) = R*x(1:2)';
        x(3:4) = R*x(3:4)';
        x(5:6) = R*x(5:6)';
        %assert((norm(x(1:2))-dold)<1e-14);
          %x(5) = 0;
        %xava=1;
    end
    
    
    % The current binding configuration state
    state = 1:p;
    
    switchtime = zeros(1, p);
    
    
    % initialize close beads bound
    for k=1:p
        for l=(k+1):p
            if state(k)==k && state(l)==l && norm(x(:,k)-x(:,l))<1e-4
                state(k)=l;
                state(l)=k;
                b = exprnd(break_mean);
                switchtime(k) = b;
                switchtime(l) = b;
            end
        end
    end
    
    v = zeros(d, p);
    
    pos_t = zeros(d, p, 10000);
    state_t = zeros(p, 10000);
    pos_t(:, :, 1) = x;
    state_t(:, 1) = state;
    
    x0 = x(:);
    
    c = 1;
    
    [rowm, colm] = ind2sub4up(1:nb);
    
    
    for step=1:nt
        if ~old_bonds
            % update bonding, handling breaking/forming within the timestep
            
            % the matrix that represents fraction of time spent bound to
            % other beads
            binding_time_fraction = zeros(p, p);
            
            % first, handle the beads that come in to timestep bound
            for k=1:p
                if state(k) ~= k
                    if switchtime(k)<(t+dt)
                        % break the bond
                        state(k)=k;
                        % register the fraction of time it spent bound
                        binding_time_fraction(k, state(k)) = (switchtime(k)-t)/dt;
                    else
                        % it stays bound the whole step
                        binding_time_fraction(k, state(k)) = 1;
                    end
                end
            end
        end
        
        
        % Update bead positions
        for k=1:p
            if state(k)==0
                partner = k;
            else
                partner = state(k);
            end
            v(:,k) = kv*(x(:,partner)-x(:,k)) - kc*x(:,k)*(norm(x(:,k)).^2);
            % add in excluded volume
            for l=1:p
                if l==k
                    continue
                end
                v(:,k) = v(:,k) + a_ev * (x(:,k) - x(:,l))...
                   .* exp(-norm(x(:,k)-x(:,l)).^2 / c_ev);
            end
        end
        
        % update variables
        drift = dt*(v/phi);
        
        sv = mod(step, 1000)+1;
        if sv == 2
            noiseArr = normrnd(0, 1, [d,p,1000]);
            expNoiseArr = exprnd(1, [nb, 1000]);
        end
        noise = sqrt(2*eps*dt)*noiseArr(:,:,sv);
        
        x = x + drift + noise;
        
        
        if fix_beads
            % rotate system to maintain symmetry fixing
            xv = x(5:6);
            r = norm(xv);
            c = xv(2)/r;
            s = xv(1)/r;
            R = [c -s; s c];
            x(1:2) = R*x(1:2)';
            x(3:4) = R*x(3:4)';
            x(5:6) = R*x(5:6)';
            %assert((norm(x(1:2))-dold)<1e-14);
              %x(5) = 0;
            %xava=1;
        end
        
        if (c > size(pos_t, 3))
            sznew = 2*size(pos_t, 3);
            if (sznew > buffer_size)
                % too long - reset trajectory
                pos_t = zeros(d, p, 10000);
                state_t = zeros(p, 10000);
                c=1;
            else
                % preallocate twice the current size
                pos_t(:,:,sznew) = 0;
                state_t(:, sznew) = 0;
            end
        end
        
        pos_t(:, :, c) = x;
        state_t(:,c) = state;
        c=c+1;
        
        t = t+dt;
        
        %D = pdist(x');
        D = zeros(1, p);
        cc = 1;
        for i=1:p
            for j=(i+1):p
                D(cc) = norm(x(:,i)-x(:,j));
                cc = cc+1;
            end
        end
        
        % check halting condition if applicable
        if ~isempty(stop)
            stopflag = 0;
            % check the stop conditions
            if planestop
                for i=1:length(stop.planes)
                    if sign(dot(x(:), stop.planes(i).normal)-stop.planes(i).bias)~= ...
                            sign(dot(x0, stop.planes(i).normal)-stop.planes(i).bias)
                        % we've crossed this plane
                        saddle_dist = norm(x(:) - stop.planes(i).v);
                        stopflag = 1;
                    end
                end
            else
                % check if beads 1 and 2 are the closest
                mindist = min(D);
                if D(1) > 1 && mindist < 0.3 
                    stopflag = 1;
                end
            end
            
            if stopflag
                pos_t = pos_t(:, :, 1:c-1);
                state_t = state_t(:, 1:c-1);
                return
            end
        end
        
        
        if old_bonds
            % Update binding state

            % form bonds
            %form_times = exprnd(eps./a(D));
            form_times = (eps./a(D)).*(expNoiseArr(:, sv)');
            [sorted, I] = sort(form_times);

            %[row, col] = ind2sub4up(I);
            row = rowm(I);
            col = colm(I);

            % 
            for j=1:length(sorted)
                % only go through times that fall in the current
                % timestep
                form_time = sorted(j);
                if form_time > dt
                    break;
                end

                b1 = row(j);
                b2 = col(j);

                assert(b1~=b2);
                if state(b1) == b1 && state(b2) == b2
                    % form a bond between these two
                    %fprintf("Forming a bond between %d and %d.\n", b1, b2);
                    state(b1) = b2;
                    state(b2) = b1;
                    % draw a breaking time for the bond
                    break_time = t + form_time + exprnd(break_mean);
                    switchtime(b1) = break_time;
                    switchtime(b2) = break_time;
                end
            end

            % break bonds
            for k=1:p
                if state(k) ~= k && switchtime(k)<t
                    state(k)=k;
                end
            end
        end
        
    end
    pos_t = pos_t(:, :, 1:c-1);
    state_t = state_t(:, 1:c-1);
end