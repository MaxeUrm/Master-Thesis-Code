%% Finer discretization error improvement plots - LW1FF

% We look at the domain \Omega = B_{1}(0) and we choose the functions

 u = @(x,t) exp(-t) * ((norm(x)^2 + x(1)*x(2)));

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

% The LW boundary conditions allow a different choice for w_Gamma:

w_gamma = @(x) exp(-t) * norm(x)^2;

% We want to see a lower error when we use the FF systems. Therefore we
% look at a system that evolves a lot different on the bdry than in the
% bulk. Therefore we have to set completly different values to the
% constants. We will generate two plots: one for finer timesteps and one
% for more dofs on the bdry.

%% Finer timesteps

% We want to look at the following values for l

l_array = [1;2;4;8;16;32];

%% Finer spartial discretizations

% We use equidistant nodes on the bdry with different dof

dof_array = [80;90;95;100;105;110;115;120;125;130;135;140;145;150;155;160;165;170;180;190;200];

%% Benchmark errors

% As a benchmark we use the AC1 System with the same constants as chosen
% above. We compare the other errors to this.

% [h_bm, tau_bm, err_L2_bm, err_H1_bm, L2_int_bm, H1_int_bm] = convergence_LW1(0.03,1000,u);

%% Compute errors

L2_l_matrix = zeros(length(l_array),3);
H1_l_matrix = zeros(length(l_array),3);
L2_l_matrix_bdry = zeros(length(l_array),3);
H1_l_matrix_bdry = zeros(length(l_array),3);

L2_dof_matrix = zeros(length(dof_array),3);
H1_dof_matrix = zeros(length(dof_array),3);
L2_dof_matrix_bdry = zeros(length(dof_array),3);
H1_dof_matrix_bdry = zeros(length(dof_array),3);

% error with rising l
%{
% dof = 200
for i = 1:length(l_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(200,l_array(i),u);
    
    L2_l_matrix(i,1) = L2_int;
    H1_l_matrix(i,1) = H1_int;
    L2_l_matrix_bdry(i,1) = L2_int_bdry;
    H1_l_matrix_bdry(i,1) = H1_int_bdry;
    
end


% dof = 400
for i = 1:length(l_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(400,l_array(i),u);

    L2_l_matrix(i,2) = L2_int;
    H1_l_matrix(i,2) = H1_int;
    L2_l_matrix_bdry(i,2) = L2_int_bdry;
    H1_l_matrix_bdry(i,2) = H1_int_bdry;
    
end

% dof = 800
for i = 1:length(l_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(800,l_array(i),u);

    L2_l_matrix(i,3) = L2_int;
    H1_l_matrix(i,3) = H1_int;
    L2_l_matrix_bdry(i,3) = L2_int_bdry;
    H1_l_matrix_bdry(i,3) = H1_int_bdry;
    
end
%}

% error with rising dof on the bdry
% l = 1
for i = 1:length(dof_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(dof_array(i),1,u);

    L2_dof_matrix(i,1) = L2_int;
    H1_dof_matrix(i,1) = H1_int;
    L2_dof_matrix_bdry(i,1) = L2_int_bdry;
    H1_dof_matrix_bdry(i,1) = H1_int_bdry;
    
end

% l = 4
for i = 1:length(dof_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(dof_array(i),4,u);

    L2_dof_matrix(i,2) = L2_int;
    H1_dof_matrix(i,2) = H1_int;
    L2_dof_matrix_bdry(i,2) = L2_int_bdry;
    H1_dof_matrix_bdry(i,2) = H1_int_bdry;
    
end

% l = 16
for i = 1:length(dof_array)
    [L2_int_bdry, H1_int_bdry, L2_int, H1_int] = convergence_LW1FF(dof_array(i),16,u);

    L2_dof_matrix(i,3) = L2_int;
    H1_dof_matrix(i,3) = H1_int;
    L2_dof_matrix_bdry(i,3) = L2_int_bdry;
    H1_dof_matrix_bdry(i,3) = H1_int_bdry;
    
end


%% Plot errors


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
subplot(3,2,1)
loglog(l_array,L2_l_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)


title('Feinere Zeitdiskretisierung','FontSize',fontSize)
xlabel('$l$','FontSize',fontSize)
%ylabel('$L^{2}$-Error','FontSize',fontSize)
%legend('$\tau = 0.1$','$\tau = 0.05$','$\tau = 0.025$','$\tau = 0.0125$','$\tau = 0.00625$','$\tau = 0.003125$','$\tau = 0.0015625$','$\tau = 0.00078125$', '$\mathcal(O)(h)$','Location','southeast','FontSize',12,'NumColumns',2)
hold off


subplot(3,2,2)
loglog(l_array,L2_l_matrix(:,2),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,2),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,2),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,2),"-x",'LineWidth',2, 'Markersize', 10)


subplot(3,2,3)
loglog(l_array,L2_l_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)


subplot(3,2,4)
loglog(l_array,L2_l_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)


subplot(3,2,5)
loglog(l_array,L2_l_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)


subplot(3,2,6)
loglog(l_array,L2_l_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
loglog(l_array,L2_l_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

loglog(l_array,H1_l_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
loglog(l_array,H1_l_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)



semilogy(dof_array,L2_dof_matrix(:,1),"-o",'LineWidth',2, 'Markersize', 10)
hold on
semilogy(dof_array,L2_dof_matrix_bdry(:,1),"-o",'LineWidth',2, 'Markersize', 10)

semilogy(dof_array,H1_dof_matrix(:,1),"-x",'LineWidth',2, 'Markersize', 10)
semilogy(dof_array,H1_dof_matrix_bdry(:,1),"-x",'LineWidth',2, 'Markersize', 10)

plot([128,128],[10^(-5),10^5],'r--')


title('Fehler bei feiner werdender Zeitdiskretisierung','FontSize',fontSize)
xlabel('DOF am Rand','FontSize',fontSize)
ylabel('Fehler','FontSize',fontSize)

legend('$L^2(\Omega)$ Fehler','$L^2(\Gamma)$ Fehler','$H^1(\Omega)$ Fehler','$H^1(\Gamma)$ Fehler','DOF = 128')

