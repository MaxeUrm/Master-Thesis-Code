%% Time_stepping_LW1 test for CH-equation with LW Boundary values

% Order 1
tic

%% Discretization bulk
% Circle
%{
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.15, [-1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

% Square
%{
fd = @(p) -min(min(min(1+p(:,2),1-p(:,2)),1+p(:,1)),1-p(:,1));
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.1, [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

% BigCircle
%
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.07, [-1,-1;1,1] );
Nodes_bulk = 5 * Nodes_bulk;
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

N_Omega = length(Nodes_bulk);
No_Elements = length(Elements_bulk);

%% Discretization surface
% Circle
%{
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

% Square
%{
Nodes_surf = zeros(1,2); % Extract boundary nodes
N_Gamma = 0;
position_array = 0;
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    if abs(x(1)) > 1 - 1/N_Gamma | abs(x(2)) > 1 - 1/N_Gamma
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

% Big Circle
%
Nodes_surf = zeros(1,2); % Extract boundary nodes
N_Gamma = 0;
position_array = 0;
for i = 1:N_Omega
    if norm(Nodes_bulk(i,:)) > 5-1/N_Omega
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

%% Starting values

% Circle
% Cone
%{
f = @(x) 2 * norm(x) -1;
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

xi_0 = ones(N_Gamma,1);

% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off
%}

% 1-function
%{
alpha_0 = ones(N_Omega,1);
xi_0 = ones(N_Gamma,1);
%}

% Peak polynom
%{
f = @(x) (norm(x)-1)^2;
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

xi_0 = zeros(N_Gamma,1);

% Plot
figure(1)
hold on
%
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end


trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off
%}

% Peak cos
%{
f = @(x) cos(pi*norm(x));
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off

xi_0 = - ones(N_Gamma,1);
%}

% Wave cos
%{
f = @(x) cos(2*pi*norm(x));
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off

xi_0 = ones(N_Gamma,1);
%}


% Square
% 1-function
%{
alpha_0 = ones(N_Omega,1);
xi_0 = ones(N_Gamma,1);
%}

% peak polynom
%{
f = @(x) ((x(1).^2 -1).^2) * ((x(2).^2 -1).^2);
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end
% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off
xi_0 = zeros(N_Gamma,1);
%}

% Peak cos
%{
f = @(x) cos(pi * x(1)) * cos(pi * x(2));
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end
% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off

xi_0 = 0;
for i = 1:N_Gamma
    xi_0 = [xi_0; f(Nodes_surf(i,:))];
end
xi_0 = xi_0(2:end,:);
%}

% Example 
%{
f = @(x) cos(4 * pi * x(1)) * cos(4 * pi * x(2));
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end
% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off

xi_0 = zeros(N_Gamma,1);
%}


% BigCircle
% Wave cos
%{
f = @(x) cos(pi*norm(x));
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = Nodes_bulk(i,:);
    alpha_0(i) = f(x);
end

% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
hold off

xi_0 = - ones(N_Gamma,1);
%}

% RIV
%
rng(42069)
alpha_0 = zeros(N_Omega,1);
for i = 1:N_Omega
    x = randn(1);
    x = 2 * (exp(5*x)/(exp(5*x)+1)) - 1; 
    alpha_0(i) = x;
end

% Plot
figure(1)
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],...
        'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off

xi_0 = zeros(N_Gamma,1);

for i = 1:N_Gamma
    xi_0(i) = alpha_0(position_array(i));
end
%}

%% Matrices
[S_bulk,M_bulk] = assembly_bulk(Nodes_bulk, Elements_bulk);
[S_surf,M_surf] = assembly_surface(Nodes_surf, Elements_surf);

%% Time discretization
T = 1;
N = 1000;

%% Coefficients
epsilon = 0.3;
sigma = 1;
delta = 0.5;
kappa = 1;

%% Compute solution
[sol_alpha,sol_beta,sol_gamma,sol_xi,sol_zeta] = time_stepping_LW1(S_bulk,M_bulk,S_surf,M_surf,Nodes_bulk, Elements_bulk, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, epsilon, sigma, delta, kappa);

%% Plotting the solution with fill3 / trisurf

% t = 0.01 

%{
figure(2)
alpha_0 = sol_alpha(:,6);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%{
figure(3)
alpha_0 = sol_alpha(:,21);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%

% t = 0.025 
%
figure(4)
alpha_0 = sol_alpha(:,51);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%

% t = 0.5
figure(5)
alpha_0 = sol_alpha(:,101);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%

%
figure(6)
alpha_0 = sol_alpha(:,151);
hold on
shading interp
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%{
figure(7)
alpha_0 = sol_alpha(:,501);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%{
figure(8)
alpha_0 = sol_alpha(:,1001);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%{
figure(9)
alpha_0 = sol_alpha(:,1501);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%{
figure(19)
alpha_0 = sol_alpha(:,2001);
hold on
%{
for i = 1:No_Elements
    fill3([Nodes_bulk(Elements_bulk(i,1),1), Nodes_bulk(Elements_bulk(i,2),1), Nodes_bulk(Elements_bulk(i,3),1)],...
        [Nodes_bulk(Elements_bulk(i,1),2), Nodes_bulk(Elements_bulk(i,2),2), Nodes_bulk(Elements_bulk(i,3),2)],...
        [alpha_0(Elements_bulk(i,1)), alpha_0(Elements_bulk(i,2)), alpha_0(Elements_bulk(i,3))],'r')
end
%}

trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none')
colorbar
shading interp
hold off
%}

%
figure(69)
subplot(1,2,1)
hold on
alpha_0 = sol_alpha(:,50);
trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none');
colorbar
shading interp
view(0,90)
hold off

subplot(1,2,2)
hold on
alpha_0 = sol_alpha(:,400);
trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none');
colorbar
shading interp
view(0,90)
hold off

%}

%% Video file creation

video = VideoWriter('Abspann','MPEG-4');
 
open(video);

for i = 1:(N+1)
    
    i
    
    alpha_0 = sol_alpha(:,i);
    trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none');
    colorbar
    shading interp
    view(0,90)
    title('Vielen Dank für Ihre Aufmerksamkeit','FontSize',24)
    frame = getframe(gcf);
    writeVideo(video,frame);
end

close(video);


%% Mass conservation
loyolagreen = 1/255*[0,104,87];
limegreen = 1/255*[132,186,91];
trueblack = 1/255*[0,0,0];
fontSize = 20;

set (0, 'defaultaxesfontname', 'Times') 
set (0, 'defaultaxesfontsize', 20)
set (0, 'defaulttextfontname', 'Times') 
set (0, 'defaulttextfontsize', 20);
set (0, 'defaultlinelinewidth', 2);
set(0, 'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

% inner mass
t = zeros(N+1,1);
m1 = zeros(N+1,1);
N = 5000;
T = 5;

for i = 1:(N+1)
    t(i) = (i-1) * (T/N); % time
    
    m1(i) = integral_approx2d(sol_alpha(:,i), Nodes_bulk, Elements_bulk);
end

figure(10)
subplot(1,2,1)
plot(t,m1,'Color',loyolagreen,'LineWidth',2)
hold on

% boundary mass
m2 = zeros(N+1,1);

for i = 1:(N+1)
    m2(i) = integral_approx1d(sol_xi(:,i), Nodes_surf, Elements_surf);
end

plot(t,m2,'Color',limegreen,'LineWidth',2)

% combined mass
plot(t,m1+m2,'Color',trueblack,'LineWidth',2)
title('Massenerhaltung','FontSize',fontSize)
xlabel('Zeit','FontSize',fontSize)
ylabel('Masse','FontSize',fontSize)
xlim([0,5])
ylim([-2,2])
legend('Masse im Inneren','Masse am Rand','Gesamtmasse','FontSize',20)
hold off

%% Energy dissipation

e = zeros(N+1,1);
for i = 1:(N+1)
    e(i) = energy_bulk(sol_alpha(:,i), epsilon, Nodes_bulk, Elements_bulk, S_bulk);
    e(i) = e(i) + energy_surf(sol_xi(:,i),delta,kappa,Nodes_surf,Elements_surf,S_surf);
end

subplot(1,2,2)
plot(t,e,'Color',loyolagreen,'LineWidth',2)
hold on
title('Energiedissipation','FontSize',fontSize)
xlabel('Zeit','FontSize',fontSize)
ylabel('Gesamtenergie','FontSize',fontSize)
xlim([0,5])
ylim([0,200])
legend('Gesamtenergie','FontSize',20)
hold off

toc
