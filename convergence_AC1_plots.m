%% Convergence plots for CH AC1


%% Solution

% We look at the domain \Omega = B_{1}(0) and we choose the functions

 u = @(x,t) exp(-t) * norm(x)^2;

% u = @(x,t) exp(-t);

% u = @(x,t) 1; 

% So we get the normal derivative

del_u = @(x,t) 2 * exp(-t);

% And u on Gamma as

u_gamma = @(x,t) exp(-t);

% We choose the chemical potential as

 w = @(x,t) exp(-t) * (norm(x)^2 -1)^2;

% w = @(x,t) exp(-t);

% w = @(x,t) 1;

% as solutions to the Cahn-Hilliard equation with HN boundary values.
% Before compution, please be sure, that the modified rhs in
% time_stepping_HN1 is in use. We now vary mesh width and timestep size and
% get back some errors in L2 and H1 norm.

%% Discretization sizes

h_vec = [0.15; 0.1; 0.07; 0.05; 0.035; 0.025]; % doubling the dof

% dofs
% 0.1 -> 362, 0.07 -> 737, 0.05 -> 1452, 0.035 -> 2955, 0.025 -> 5809,
% 0.019 -> 10045, 0.013 -> 21465

tau_vec = [10; 20; 40; 80; 160; 320; 640; 1280]; % doubling the number of timesteps

h_arr = zeros(5,1);
tau_arr = zeros(8,1);

L2_err_matrix = zeros(5,8);
H1_err_matrix = zeros(5,8);

%% Compute errors
%

for i = 1:6
    for j = 1:8
        [h, tau, err_L2, err_H1, L2_int, H1_int] = convergence_AC1(h_vec(i),tau_vec(j),u);
        
        i
        j
        
        L2_err_matrix(i,j) = sqrt(L2_int);
        H1_err_matrix(i,j) = sqrt(H1_int);
        h_arr(i) = h;
        tau_arr(j) = tau;
    end
end
%}

%% Plotting space

loyolagreen = 1/255*[0,104,87];
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

figure(1) 
subplot(2,2,1)
loglog(h_arr,L2_err_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(h_arr,L2_err_matrix(:,2),"-*",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,3),"-+",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,4),"-x",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,5),"-v",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,6),"-diamond",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,7),"-square",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,L2_err_matrix(:,8),"-^",'LineWidth',2, 'Markersize', 10)

temp = 0.05 * h_arr;

loglog(h_arr,temp,"--",'LineWidth',2)

title('$L^{2}$ Fehler','FontSize',fontSize)
xlabel('$h$','FontSize',fontSize)
ylabel('$L^{2}$-Error','FontSize',fontSize)
ylim([10^(-4),1])
legend('$\tau = 0.1$','$\tau = 0.05$','$\tau = 0.025$','$\tau = 0.0125$','$\tau = 0.00625$','$\tau = 0.003125$','$\tau = 0.0015625$','$\tau = 0.00078125$', '$\mathcal(O)(h)$','Location','southeast','FontSize',12,'NumColumns',2)
hold off

subplot(2,2,2)
loglog(h_arr,H1_err_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(h_arr,H1_err_matrix(:,2),"-*",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,3),"-+",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,4),"-x",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,5),"-v",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,6),"-diamond",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,7),"-square",'LineWidth',2, 'Markersize', 10)
loglog(h_arr,H1_err_matrix(:,8),"-^",'LineWidth',2, 'Markersize', 10)

temp = 0.17 * h_arr;

loglog(h_arr,temp,"--",'LineWidth',2)

title('$H^{1}$ Fehler','FontSize',fontSize)
xlabel('$h$','FontSize',fontSize)
ylabel('$H^{1}$-Error','FontSize',fontSize)
ylim([10^(-4),1])
legend('$\tau = 0.1$','$\tau = 0.05$','$\tau = 0.025$','$\tau = 0.0125$','$\tau = 0.00625$','$\tau = 0.003125$','$\tau = 0.0015625$','$\tau = 0.00078125$', '$\mathcal(O)(h)$','Location','southeast','FontSize',12,'NumColumns',2)

hold off

%% Plotting time

loyolagreen = 1/255*[0,104,87];
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

subplot(2,2,3)
loglog(tau_arr,L2_err_matrix(1,:),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(tau_arr,L2_err_matrix(2,:),"-*",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,L2_err_matrix(3,:),"-+",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,L2_err_matrix(4,:),"-x",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,L2_err_matrix(5,:),"-diamond",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,L2_err_matrix(6,:),"-square",'LineWidth',2, 'Markersize', 10)

temp = 0.7 * tau_arr;

loglog(tau_arr,temp,"--",'LineWidth',2)

title('$L^{2}$ Fehler','FontSize',fontSize)
xlabel('$\tau$','FontSize',fontSize)
ylabel('$L^{2}$-Error','FontSize',fontSize)
ylim([10^(-4),1])
legend('dof 131','dof 362','dof 737','dof 1452','dof 2955','dof 5809', '$\mathcal(O)(\tau)$','Location','southeast','FontSize',12,'NumColumns',2)
hold off

subplot(2,2,4)
loglog(tau_arr,H1_err_matrix(1,:),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(tau_arr,H1_err_matrix(2,:),"-*",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,H1_err_matrix(3,:),"-+",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,H1_err_matrix(4,:),"-x",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,H1_err_matrix(5,:),"-diamond",'LineWidth',2, 'Markersize', 10)
loglog(tau_arr,H1_err_matrix(6,:),"-square",'LineWidth',2, 'Markersize', 10)

temp = 2.8 * tau_arr;

loglog(tau_arr,temp,"--",'LineWidth',2)

title('$H^{1}$ Fehler','FontSize',fontSize)
xlabel('$\tau$','FontSize',fontSize)
ylabel('$H^{1}$-Error','FontSize',fontSize)
ylim([10^(-4),1])
legend('dof 131','dof 362','dof 737','dof 1452','dof 2955','dof 5809', '$\mathcal(O)(\tau)$','Location','southeast','FontSize',12,'NumColumns',2)
hold off
