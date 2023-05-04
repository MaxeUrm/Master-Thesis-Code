%% Convergence plot function
% This function computes a solution to the manipulated CH-System which can
% be comared to the true solution

function [h, tau, err_L2, err_H1, L2_int, H1_int] = convergence_LW2(mesh,N,u,nr)

%% Discretization bulk

% Choose meshsize
%
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, mesh, [-1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

N_Omega = length(Nodes_bulk);
No_Elements = length(Elements_bulk);

h = max_width(Nodes_bulk,Elements_bulk,No_Elements);

%% Discretization surface
% Circle
%
Nodes_surf = zeros(1,2); % Extract boundary nodes
N_Gamma = 0;
position_array = 0;
for i = 1:N_Omega
    if norm(Nodes_bulk(i,:)) > 1-1/N_Omega
        Nodes_surf = [Nodes_surf; Nodes_bulk(i,:)];
        N_Gamma = N_Gamma + 1; % count for N_Gamma
        position_array = [position_array; i];
        
    end
end
Nodes_surf = Nodes_surf(2:(N_Gamma+1),:);
position_array = position_array(2:(N_Gamma+1),:); % for extracting boundary elements

Elements_surf = zeros(1,2); % Extract boundary nodes
for i = 1:(length(Elements_bulk))
    x1 = Elements_bulk(i,1); % Element i, Node 1
    x2 = Elements_bulk(i,2); % Element i, Node 2
    x3 = Elements_bulk(i,3); % Element i, Node 3
    
    count = 0;
    vec = 0;
    
    if ismember(x1,position_array)
        [a,b] = ismember(x1,position_array); % b is new place in Nodes_surf
        count = count +1;
        vec = [vec,b];
    end
    
    if ismember(x2,position_array)
        [a,b] = ismember(x2,position_array);
        count = count +1;
        vec = [vec,b];
    end
    
    if ismember(x3,position_array)
        [a,b] = ismember(x3,position_array);
        count = count +1;
        vec = [vec,b];
    end
    
    vec = vec(2:end);
    
    
    if count == 2 % two nodes lie on the boundary
        Elements_surf = [Elements_surf; vec];
    end
end

Elements_surf = Elements_surf(2:end,:);
%}

%% Starting value 

alpha_0 = ones(N_Omega,1);

% Peak norm
%
f = @(x) norm(x)^2;
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

% Plot
%{
figure(1)
hold on
trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
%shading interp
hold off
%}

xi_0 = ones(N_Gamma,1);

%% Matrices
[S_bulk,M_bulk] = assembly_bulk(Nodes_bulk, Elements_bulk);
[S_surf,M_surf] = assembly_surface(Nodes_surf, Elements_surf);

%% Time discretization
T = 1;
tau = T/N;

%% Coefficients
epsilon = 1;
sigma = 1;
delta = 1;
kappa = 1;

%% Compute solution
[sol_alpha,sol_beta,sol_gamma,sol_xi,sol_zeta] = time_stepping_LW2(S_bulk,M_bulk,S_surf,M_surf,Nodes_bulk, Elements_bulk, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, epsilon, sigma, delta, kappa, nr);

%% Error estimate

% true solutions

u_sol = alpha_0;

for i = 1:N
    v = zeros(N_Omega,1);
    for j = 1:N_Omega
        v(j) = u(Nodes_bulk(j,:),i*tau);
    end
    u_sol = [u_sol,v];
end

u_sol = u_sol - sol_alpha; % Difference on every timestep

err_L2 = zeros(N+1,1);
L2_int = 0;

err_H1 = zeros(N+1,1);
H1_int = 0;

for i = 2:(N+1)
    err_L2(i) = transpose(u_sol(:,i)) * M_bulk * u_sol(:,i);
    L2_int = L2_int + (err_L2(i)+err_L2(i-1)) * tau/2;
    
    err_H1(i) = transpose(u_sol(:,i)) * S_bulk * u_sol(:,i);
    H1_int = H1_int + (err_H1(i)+err_H1(i-1)) * tau/2;
end

%{
figure(2)
alpha_0 = sol_alpha(:,end);
hold on

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}
end