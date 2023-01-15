%% Convergence plot function
% This function computes a solution to the manipulated CH-System which can
% be comared to the true solution

function [h, tau, err_L2, err_H1, L2_int, H1_int] = convergence_HN1(mesh,N,u)

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

%% Starting value 

alpha_0 = ones(N_Omega,1);

% Peak norm
%
f = @(x) (norm(x)^2 -1)^2;
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



%% Matrices
[S_bulk,M_bulk] = assembly_bulk(Nodes_bulk, Elements_bulk);

%% Time discretization
T = 1;
tau = T/N;

%% Coefficients
epsilon = 0.3;
sigma = 1;

%% Compute solution with manipulated rhs
[sol_alpha,sol_beta] = time_stepping_HN1(S_bulk,M_bulk,Nodes_bulk, Elements_bulk, alpha_0, T, N, epsilon, sigma);

%% Video file creation
%{
video = VideoWriter('test','MPEG-4');

open(video);

for i = 1:(N+1)
    
    %i
    
    alpha_0 = sol_alpha(:,i);
    trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none');
    colorbar
    shading interp
    view(0,90)
    frame = getframe(gcf);
    writeVideo(video,frame);
end

close(video);
%}

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
    L2_int = L2_int + err_L2(i) * tau;
    
    err_H1(i) = transpose(u_sol(:,i)) * S_bulk * u_sol(:,i);
    H1_int = H1_int + err_H1(i) * tau;
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