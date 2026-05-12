%% 2D CHROMATIN POLYMER WITH CROSSLINKING
clear; clc; close all;

%% ---------------- PARAMETERS ----------------

N = 40;                     % beads
steps = 50000;
dt = 1e-3;

zeta = 2.5e-3;
kBT  = 4.1;

Lp = 50;
Nk = 17;
R0 = Nk*(2*Lp);

alpha = 0.2176;

% excluded volume
cEV = 8.305e-5;
aEV = 3.268e-5;

% nuclear radius
R_nuc = 500;

% crosslink parameters
bind_dist = 90;            % nm
k_cross = 10;              % crosslink spring
p_on  = 0.001;             % binding probability
p_off = 0.0001;            % unbinding probability

noise = sqrt(2*kBT*dt/zeta);

%% ---------------- INITIAL POLYMER ----------------

theta = linspace(0,pi,N);

X = [R_nuc*cos(theta)', R_nuc*sin(theta)'];

% tether ends to boundary
X(1,:) = [-R_nuc,0];
X(N,:) = [R_nuc,0];

crosslinks = zeros(N);     % adjacency matrix

%% ---------------- PLOT SETUP ----------------

figure('color','w')
hold on

th = linspace(0,2*pi,200);
plot(R_nuc*cos(th),R_nuc*sin(th),'k')

chain = plot(X(:,1),X(:,2),'o-','LineWidth',2);

axis equal
xlim([-R_nuc R_nuc])
ylim([-R_nuc R_nuc])

title('Chromatin Polymer Simulation')

%% ---------------- SIMULATION ----------------

for t = 1:steps

F = zeros(N,2);

%% ---------- WLC NEIGHBOR FORCES ----------

for i = 1:N-1

    rij = X(i+1,:) - X(i,:);
    r = norm(rij);

    FWLC = alpha * (-1 + 1/(1-r/R0)^2 + 4*r/R0);

    fvec = FWLC*(rij/r);

    F(i,:)   = F(i,:) + fvec;
    F(i+1,:) = F(i+1,:) - fvec;

end

%% ---------- EXCLUDED VOLUME ----------

for i = 1:N
for j = i+1:N

    rij = X(i,:) - X(j,:);
    r2 = sum(rij.^2);

    f = cEV * rij * exp(-aEV*r2);

    F(i,:) = F(i,:) + f;
    F(j,:) = F(j,:) - f;

end
end

%% ---------- CROSSLINK FORCES ----------

for i = 1:N
for j = i+2:N

    rij = X(i,:) - X(j,:);
    r = norm(rij);

    % binding
    if crosslinks(i,j)==0 && r < bind_dist
        if rand < p_on
            crosslinks(i,j) = 1;
            crosslinks(j,i) = 1;
        end
    end

    % unbinding
    if crosslinks(i,j)==1
        if rand < p_off
            crosslinks(i,j) = 0;
            crosslinks(j,i) = 0;
        else

            f = -k_cross*(r-30)*(rij/r);

            F(i,:) = F(i,:) + f;
            F(j,:) = F(j,:) - f;

        end
    end

end
end

%% ---------- BROWNIAN FORCE ----------

xi = noise * randn(N,2);

%% ---------- UPDATE ----------

X = X + (dt/zeta)*F + xi;

%% ---------- NUCLEAR BOUNDARY ----------

for i = 1:N

    r = norm(X(i,:));

    if r > R_nuc
        X(i,:) = X(i,:) * (R_nuc/r);
    end

end

%% ---------- TETHERS ----------

X(1,:) = [-R_nuc,0];
X(N,:) = [R_nuc,0];

%% ---------- PLOT ----------

if mod(t,50)==0
    set(chain,'XData',X(:,1),'YData',X(:,2))
    drawnow
end

end