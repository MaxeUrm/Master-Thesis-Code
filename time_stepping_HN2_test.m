%% time_stepping_HN2 test for CH-equation with HN boundary values

% Order 2

%% Discretization bulk

set (0, 'defaultaxesfontname', 'Times') 
set (0, 'defaultaxesfontsize', 20)
set (0, 'defaulttextfontname', 'Times') 
set (0, 'defaulttextfontsize', 20);
set(0, 'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

% Circle
%{
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.025, [-1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

% Square
%{
fd = @(p) -min(min(min(1+p(:,2),1-p(:,2)),1+p(:,1)),1-p(:,1));
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.2, [-1,-1;1,1], [-1,-1;-1,1;1,-1;1,1] );
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}


% BigCircle
%
fd = @(p) sqrt(sum(p.^2,2)) - 1;
fh = @(p) ones(size(p,1),1);
[Nodes_bulk,Elements_bulk, stats] = distmesh( fd, fh, 0.025, [-1,-1;1,1] );
Nodes_bulk = 5 * Nodes_bulk;
% patch( 'vertices', Nodes_bulk, 'faces', Elements_bulk, 'facecolor', [.9, .9, .9] )
%}

N_Omega = length(Nodes_bulk);
No_Elements = length(Elements_bulk);

h = max_width(Nodes_bulk,Elements_bulk,No_Elements)

%% Starting values

% 1-function
%
alpha_0 = ones(N_Omega,1);
%

% Circle

% Peak polynom
%{
f = @(x) (norm(x)-1)^2;
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
hold off
%}

% Peak cos
%{
f = @(x) cos(pi*norm(x));
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
hold off
%}
trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0)
colorbar
hold off
%}

% Wave cos
%{
f = @(x) cos(2*pi*norm(x));
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
hold off
%}


% Square

% 0-function
%{
alpha_0 = zeros(N_Omega,1);
%}

% peak polynom
%{
f = @(x) ((x(1).^2 -1).^2) * ((x(2).^2 -1).^2);
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
%}

% Peak cos
%{
f = @(x) cos(pi * x(1)) * cos(pi * x(2));
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
%}

% Example 
%{
f = @(x) cos(4 * pi * x(1)) * cos(4 * pi * x(2));
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
hold off
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
colormap('jet')
pbaspect([2 2 1])
colorbar
shading interp
hold off

%}


%% Matrices
[S_bulk,M_bulk] = assembly_bulk(Nodes_bulk, Elements_bulk);

%% Time discretization
T = 5;
N = 5000;

%% Coefficients
epsilon = 0.3;

sigma = 1;

%% Number of weights for Newton-Cotes formula
nr = 5;

%% Compute solution
[sol_alpha,sol_beta] = time_stepping_HN2(S_bulk,M_bulk,Nodes_bulk, Elements_bulk, alpha_0, T, N, epsilon, sigma, nr);

% solution vectors
% alpha_sol = sol(1:((N+1)*N_Omega),1);
% beta_sol = sol(((N+1)*N_Omega+1):((2*N+1)*N_Omega));

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

%{
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

video = VideoWriter('CH_HN2_vid','MPEG-4');

open(video);

for i = 1:(N+1)
    
    i
    
    alpha_0 = sol_alpha(:,i);
    trisurf(Elements_bulk, Nodes_bulk(:,1), Nodes_bulk(:,2), alpha_0, 'EdgeColor', 'none');
    colormap('jet')
    colorbar
    caxis([-1,1])
    pbaspect([2 2 1])
    shading interp
    view(0,90)
    frame = getframe(gcf);
    writeVideo(video,frame);
end

close(video);


%% Mass conservation
loyolagreen = 1/255*[0,104,87];
limegreen = 1/255*[132,186,91];
deepcarrotorange = [0.9100, 0.4100, 0.1700];
jeansblue = [0.1, 0.2, 0.6];
fontSize = 20;

set (0, 'defaultaxesfontname', 'Times') 
set (0, 'defaultaxesfontsize', 20)
set (0, 'defaulttextfontname', 'Times') 
set (0, 'defaulttextfontsize', 20);
set (0, 'defaultlinelinewidth', 2);
set(0, 'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

t = zeros(N+1,1);
m = zeros(N+1,1);

for i = 1:(N+1)
    t(i) = (i-1) * (T/N); % time
    
    m(i) = integral_approx2d(sol_alpha(:,i), Nodes_bulk, Elements_bulk);
end

figure(6)
subplot(1,2,1)
plot(t,m,'Color', loyolagreen, 'LineWidth',2)
hold on
title('Massenerhaltung','FontSize',fontSize)
xlabel('Zeit','FontSize',fontSize)
ylabel('Masse','FontSize',fontSize)
xlim([0,5])
ylim([-1,5])
legend('Masse im Inneren','FontSize',20)
hold off

%% Energy dissipation

en = zeros(N+1,1);
for i = 1:(N+1)
    en(i) = energy_bulk(sol_alpha(:,i), epsilon, Nodes_bulk, Elements_bulk, S_bulk);
end

subplot(1,2,2)
plot(t,en,'Color',loyolagreen,'LineWidth',2)
hold on
title('Energiedissipation','FontSize',fontSize)
xlabel('Zeit','FontSize',fontSize)
ylabel('Gesamtenergie','FontSize',fontSize)
xlim([0,5])
ylim([0,200])
legend('Energie im Inneren','FontSize',20)
hold off

%% Energy difference inequality
ed = zeros(N-2,1);
for i = 3:(N+1)
    ed(i) = transpose(sol_alpha(:,i)-sol_alpha(:,i-1)) * M_bulk * (sol_alpha(:,i)-sol_alpha(:,i-1)); %+ transpose(sol_xi(:,i)-sol_xi(:,i-1)) * M_surf * (sol_xi(:,i)-sol_xi(:,i-1));
end

figure(7)
semilogy(t,ed,'Color',loyolagreen,'LineWidth',2)
title('Störterme','FontSize',fontSize)
xlabel('Zeit','FontSize',fontSize)
ylabel('Energie','FontSize',fontSize)
pbaspect([3,1,1])
xlim([0,5])
legend('Störterme','FontSize',20)

