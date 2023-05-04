%% error improvement plot functions

function [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_AC1FF(N_Gamma,k,u)

%% Discretization bulk
    
% Circle
%
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.05, [-1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

N_Omega = length(Nodes_bulk);
No_Elements = length(Elements_bulk);

%% Extract boundary nodes
% Circle
%
Nodes_bdry = zeros(1,2); % Extract boundary nodes
N_bdry = 0;
position_array = 0;
for i = 1:N_Omega
    if norm(Nodes_bulk(i,:)) > 1-1/N_Omega
        Nodes_bdry = [Nodes_bdry; Nodes_bulk(i,:)];
        N_bdry = N_bdry + 1; % count for N_Gamma
        position_array = [position_array; i];
        
    end
end
Nodes_bdry = Nodes_bdry(2:(N_bdry+1),:);
position_array = position_array(2:(N_bdry+1),:); % for extracting boundary elements

Elements_bdry = zeros(1,2); % Extract boundary nodes
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
        Elements_bdry = [Elements_bdry; vec];
    end
end

Elements_bdry = Elements_bdry(2:end,:);
%}

%% Boundary discretization

% Equidistant
%{
v = linspace(0,2*pi,N_Gamma+1);
v = transpose(v(1:N_Gamma));

Nodes_surf = [cos(v),sin(v)];

Elements_surf = zeros(N_Gamma,2);
for i = 1:N_Gamma
    Elements_surf(i,:) = [i,i+1];
    if i == N_Gamma
        Elements_surf(i,:) = [N_Gamma,1];
    end
end
%}

% Identical
%
Nodes_surf = Nodes_bdry;
Elements_surf = Elements_bdry;
N_Gamma = N_bdry;
%}

[B,Bdry] = assembly_tracematrix(Nodes_bulk,Nodes_bdry,Nodes_surf,N_Omega,N_bdry,N_Gamma,position_array);

%% Starting values

%
 f = @(x) norm(x)^2 + x(1)*x(2) ;
% f = @(x) 1;
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

xi_0 = ones(N_Gamma,1);

% Plot
figure(2)
hold on

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
shading interp
colorbar
hold off
%}


 B_I = assembly_interpolmatrix(Bdry, Nodes_surf, N_Omega, N_bdry, N_Gamma);

 xi_0 = B_I * alpha_0;

%f = @(x) 1 + x(1)*x(2) ;
%{
f = @(x) 1;
for i = 1:N_Gamma
    x = Nodes_surf(i,:);
    xi_0(i) = f(x);
end
%}

%% Matrices
[S_bulk,M_bulk] = assembly_bulk(Nodes_bulk, Elements_bulk);
[S_surf,M_surf] = assembly_surface(Nodes_surf, Elements_surf);

%% Time discretization
T = 1;
N = 10;
l = 1; % refinement in the bulk

tau = T/N;

%% Coefficients
epsilon = 1; % big interaction length
sigma = 10; % slow separation
delta = 1; % small interaction length
kappa = 0.25; % fast separation

%% Compute solution
[sol_alpha,sol_beta,sol_gamma,sol_xi] = time_stepping_AC1FF(S_bulk,M_bulk,Nodes_bulk,Elements_bulk,S_surf,M_surf,B,B_I, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, k, l, epsilon, sigma, delta, kappa);

%% Error estimate

% true solutions

% bulk
u_sol = alpha_0;

for i = 1:(k*N)
    v = zeros(N_Omega,1);
    for j = 1:N_Omega
        v(j) = u(Nodes_bulk(j,:),i*tau/k);
    end
    u_sol = [u_sol,v];
end

u_sol = u_sol - sol_alpha; % Difference on every timestep

err_L2 = zeros(k*N+1,1);
L2_int = 0;

err_H1 = zeros(k*N+1,1);
H1_int = 0;

for i = 2:(k*N+1)
    err_L2(i) = transpose(u_sol(:,i)) * M_bulk * u_sol(:,i);
    L2_int = L2_int + err_L2(i) * tau/k;
    
    err_H1(i) = transpose(u_sol(:,i)) * S_bulk * u_sol(:,i);
    H1_int = H1_int + err_H1(i) * tau/k;
end


% surface
u_sol = xi_0;

for i = 1:(l*N)
    v = zeros(N_Gamma,1);
    for j = 1:N_Gamma
        v(j) = u(Nodes_surf(j,:),i*(tau/l));
    end
    u_sol = [u_sol,v];
end

u_sol = u_sol - sol_xi;

err_L2 = zeros(l*N+1,1);
L2_int_bdry = 0;

err_H1 = zeros(l*N+1,1);
H1_int_bdry = 0;

for i = 2:(l*N+1)
    err_L2(i) = transpose(u_sol(:,i)) * M_surf * u_sol(:,i);
    L2_int_bdry = L2_int_bdry + err_L2(i) * (tau/l);
    
    err_H1(i) = transpose(u_sol(:,i)) * S_surf * u_sol(:,i);
    H1_int_bdry = H1_int_bdry + err_H1(i) * (tau/l);
end
% Plot
%{
figure(10)
hold on

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), sol_alpha(:,500), 'EdgeColor', 'none')
shading interp
colorbar
hold off

figure(11)
hold on

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), sol_alpha(:,1000), 'EdgeColor', 'none')
shading interp
colorbar
hold off
%}
end