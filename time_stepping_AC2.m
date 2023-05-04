%% Time-stepping scheme slover for CH-System with Allen-Cahn-Type Boundary

% Order 2

function [sol_alpha,sol_beta,sol_gamma,sol_xi] = time_stepping_AC2(S_bulk,M_bulk,S_surf,M_surf,Nodes_bulk, Elements_bulk, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, epsilon, sigma, delta, kappa, nr)
    
    %% dof in the bulk
    N_Omega = length(Nodes_bulk);
    % dof on the surface
    N_Gamma = length(Nodes_surf);
    
    %% tau
    tau = T/N;
    
    %% assembly of the time-stepping matrix
    K = spalloc(2*(N_Omega + N_Gamma), 2*(N_Omega + N_Gamma), 10*(N_Omega + N_Gamma));
    
    B = assembly_B(Nodes_bulk, Nodes_surf);
    
        % first equation
        K(1:N_Omega,1:N_Omega) = M_bulk;
        K(1:N_Omega,(N_Omega+1):(2 * N_Omega)) = tau * sigma * S_bulk;
        
        % second equation
        K((N_Omega + 1):(2 * N_Omega),1:N_Omega) = epsilon * S_bulk;
        K((N_Omega + 1):(2 * N_Omega),(N_Omega + 1):(2 * N_Omega)) = -M_bulk;
        K((N_Omega + 1):(2 * N_Omega),(2 * N_Omega +1):(2 * N_Omega + N_Gamma)) = -epsilon * transpose(B) * M_surf;
        
        % thrid equation
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + 1):(2 * N_Omega + N_Gamma)) = epsilon * tau * M_surf;
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = M_surf + tau * delta * kappa * S_surf;
        
        % fourth equation
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),1:N_Omega) = B;
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = - diag(ones(N_Gamma,1));
    
    [ML_bulk,ML_surf] = assembly_ML(Nodes_bulk,Elements_bulk,Nodes_surf,Elements_surf,N_Omega,N_Gamma);
   
    
    % Anfangswerte
    sol_alpha = alpha_0;
    sol_beta = zeros(N_Omega,1);
    sol_gamma = zeros(N_Gamma,1);
    sol_xi = xi_0;
    
    sol_n = ones(2*(N_Omega+N_Gamma),1);
    
    % get u_1 by using the method of order 1
    %{
    alpha_n = sol_alpha(:,1);
    xi_n = sol_xi(:,1);
        
    v_1 = - M_bulk * alpha_n;
    v_2 = (1/epsilon) * ML_bulk * W_minus_bulk(alpha_n);
    v_3 = - M_surf * xi_n + (tau/delta) * ML_surf * W_minus_surf(xi_n);
    v_4 = zeros(N_Gamma,1);
        
    b = [v_1; v_2; v_3; v_4]; % rhs
    
    % true solution for CH AC1
    %{
    fun = @(x) K * x + b + f_nonlin_AC1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta); % final funtion to solve
    %}
    
    % Inhomogenous rhs for convergence plot tests OG
        %
        i = 1;
        rhs1 = @(x) exp(i*tau) * norm(x)^2 - 8 * sigma * exp(i*tau) * (2*norm(x)^2 -1);
        v1 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v1(j) = rhs1(Nodes_bulk(j,:));
        end
        
        v1 = tau * ML_bulk * v1; % Mass lumping for approximation
        
        rhs2 = @(x) - 4.0 * epsilon * exp(i*tau) + 1/epsilon * ((exp(i*tau) * norm(x)^2)^3 - exp(i*tau) * norm(x)^2) -exp(i*tau) * (norm(x)^2 -1)^2;
        v2 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v2(j) = rhs2(Nodes_bulk(j,:));
        end
        
        v2 = ML_bulk * v2; % Mass lumping for approximation
        
        rhs3 = @(x) exp(i*tau) - 4.0 * delta * kappa * exp(i*tau) * 0 + 1/delta * (exp(3*i*tau) - exp(i*tau)) + 2 * epsilon * exp(i*tau);
        v3 = zeros(N_Gamma,1);
        
        for j = 1:N_Gamma
            v3(j) = rhs3(Nodes_surf(j,:));
        end
        v3 = tau * ML_surf * v3;
        
        v = [v1;v2;v3;zeros(N_Gamma,1)];
        fun = @(x) (K * x) + b + f_nonlin_AC1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta) - v ; % final funtion to solve
        %
        
    % sol_n = fsolve(fun, sol_n)
        
    tol = 0.000001;
    sol_n = newton_solver_AC1(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, tau, delta, sol_n, tol);
    %}
        
    % update solution vector
    %{
    % alpha
    sol_alpha = [sol_alpha,sol_n(1:N_Omega)];
    % beta
    sol_beta = [sol_beta,sol_n((N_Omega + 1):(2 * N_Omega))];
    %}
    % gamma
    %sol_gamma = [sol_gamma,sol_n((2 * N_Omega + 1):(2 * N_Omega + N_Gamma))];
    %{
    % xi
    sol_xi = [sol_xi,sol_n((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma))];
    %}
    
    % manipulate the solution vectors
    %
    i=1;
    temp = zeros(N_Omega,1);
    fctn = @(x,t) exp(-10*t) * norm(x)^2;
    for i = 1:N_Omega
        temp(i) = fctn(Nodes_bulk(i,:),tau);
    end
    sol_alpha = [sol_alpha,temp];
    
    temp = zeros(N_Omega,1);
    fctn = @(x,t) exp(-10*t) * (norm(x)^2 -1)^2;
    for i = 1:N_Omega
        temp(i) = fctn(Nodes_bulk(i,:),tau);
    end
    sol_beta = [sol_beta,temp];
    
    temp = zeros(N_Gamma,1);
    fctn = @(x,t) exp(-10*t) * norm(x)^2;
    for i = 1:N_Gamma
        temp(i) = fctn(Nodes_surf(i,:),tau);
    end
    sol_xi = [sol_xi,temp];
    
    temp = zeros(N_Gamma,1);
    fctn = @(x,t) 2 * exp(-10*t);
    for i = 1:N_Gamma
        temp(i) = fctn(Nodes_surf(i,:),tau);
    end
    sol_gamma = [sol_gamma,temp];
    %}
    
    %% now order 2
    
    weights = getweights(0,1,nr);
    
    %% assembly of the time-stepping matrix
    
    K = spalloc(2*(N_Omega + N_Gamma), 2*(N_Omega + N_Gamma), 10*(N_Omega + N_Gamma));
    
    B = assembly_B(Nodes_bulk, Nodes_surf);
    
        % first equation
        K(1:N_Omega,1:N_Omega) = M_bulk;
        K(1:N_Omega,(N_Omega+1):(2 * N_Omega)) = tau * sigma * 0.5 * S_bulk;
        
        % second equation
        K((N_Omega + 1):(2 * N_Omega),1:N_Omega) = epsilon * 0.5 * S_bulk;
        K((N_Omega + 1):(2 * N_Omega),(N_Omega + 1):(2 * N_Omega)) = - 0.5 * M_bulk;
        K((N_Omega + 1):(2 * N_Omega),(2 * N_Omega +1):(2 * N_Omega + N_Gamma)) = - epsilon * 0.5 * transpose(B) * M_surf;
        
        % thrid equation
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + 1):(2 * N_Omega + N_Gamma)) = epsilon * tau * 0.5 * M_surf;
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = M_surf + tau * delta * kappa * 0.5 * S_surf;
        
        % fourth equation
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),1:N_Omega) = B;
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = - diag(ones(N_Gamma,1));
    
    
    %% loop
    for i = 2:N
        
        i % Step
        
        % construct rhs out of prevous step
        
        alpha_n = sol_alpha(:,i);
        alpha_pn = sol_alpha(:,i-1);
        
        beta_n = sol_beta(:,i);
        
        xi_n = sol_xi(:,i);
        xi_pn = sol_xi(:,i-1);
        
        gamma_n = sol_gamma(:,i);
        
        
        v_1 = - M_bulk * alpha_n + tau * sigma * 0.5 * S_bulk * beta_n;
        v_2 = epsilon * 0.5 * S_bulk * alpha_n + (1/epsilon) * (1.5 * M_bulk * W_minus_bulk(alpha_n) - 0.5 * M_bulk * W_minus_bulk(alpha_pn)) - 0.5 * M_bulk * beta_n - epsilon * 0.5 * transpose(B) * M_surf * gamma_n;
        v_3 = - M_surf * xi_n + (tau * delta * kappa * 0.5) * S_surf * xi_n + tau * (1/delta) * (1.5 * M_surf * W_minus_surf(xi_n) - 0.5 * M_surf * W_minus_surf(xi_pn)) + epsilon * tau * 0.5 * M_surf * gamma_n;
        v_4 = zeros(N_Gamma,1);
        
        b = [v_1; v_2; v_3; v_4]; % rhs
        
        % true solution for CH AC1
        %{
        fun = @(x) K * x + b + f_nonlin_AC2(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta, weights, alpha_n, xi_n); % final funtion to solve
        %}
        
        % Inhomogenous rhs for convergence plot tests OG
        %
        rhs1 = @(x,t) -10*exp(-10*t) * norm(x)^2 - 8.0 * sigma * exp(-10*t) * (2.0*norm(x)^2 -1);
        %rhs1 = @(x,t) exp(t) * cos(pi*norm(x)^2) + sigma * 4*pi*exp(t) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2));
        v1 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v1(j) = rhs1(Nodes_bulk(j,:), (i-0.5)*tau);
        end
        
        v1 = tau * M_bulk * v1; % Mass lumping for approximation
        
        rhs2 = @(x,t) - 4.0 * epsilon * exp(-10*t) + 1/epsilon * ((exp(-10*t) * norm(x)^2)^3 - exp(-10*t) * norm(x)^2) -exp(-10*t) * (norm(x)^2 -1)^2;
        %rhs2 = @(x,t) epsilon * 4*pi*exp(t) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2)) - exp(t) * cos(pi*norm(x)^2) + 1/epsilon * ((exp(t) * cos(pi*norm(x)^2))^3 - exp(t) * cos(pi*norm(x)^2));
        v2 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v2(j) = rhs2(Nodes_bulk(j,:), (i-0.5)*tau);
        end
        
        v2 = M_bulk * v2; % Mass lumping for approximation
        
        rhs3 = @(x,t) -10*exp(-10*t) + (1/delta) * ((exp(-10*t))^3 - exp(-10*t)) + 2.0 * epsilon * exp(-10*t);
        %rhs3 = @(x,t) exp(t) + 1/delta * ((exp(t) * cos(pi))^3 - exp(t) * cos(pi));
        v3 = zeros(N_Gamma,1);
        
        for j = 1:N_Gamma
            v3(j) = rhs3(Nodes_surf(j,:), (i-0.5)*tau);
        end
        v3 = tau * M_surf * v3;
        
        v = [v1;v2;v3;zeros(N_Gamma,1)];
        fun = @(x) (K * x) + b - v + f_nonlin_AC2(M_bulk, M_surf, x, N_Omega, N_Gamma, epsilon, tau, delta, weights, alpha_n, xi_n); % final funtion to solve
        %}
        
        % sol_n = fsolve(fun, sol_n)
        
        tol = 0.0000000001;
        sol_n = newton_solver_AC2(M_bulk, M_surf, K, fun, N_Omega, N_Gamma, epsilon, tau, delta, sol_n, tol, weights, alpha_n, xi_n);
        
        
        % update solution vector
        % alpha
        sol_alpha = [sol_alpha,sol_n(1:N_Omega)];
        % beta
        sol_beta = [sol_beta,sol_n((N_Omega + 1):(2 * N_Omega))];
        % gamma
        sol_gamma = [sol_gamma,sol_n((2 * N_Omega + 1):(2 * N_Omega + N_Gamma))];
        % xi
        sol_xi = [sol_xi,sol_n((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma))];
        
    end
end