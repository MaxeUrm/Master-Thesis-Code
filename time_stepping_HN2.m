%% Time-stepping scheme slover for CH-System with homogenous Neumann boundary values

% Order 2

function [sol_alpha,sol_beta] = time_stepping_HN2(S_bulk,M_bulk,Nodes_bulk, Elements_bulk, alpha_0, T, N, epsilon, sigma, nr);

    %% dof in the bulk
    N_Omega = length(Nodes_bulk);
    
    %% tau
    tau = T/N;
    
    %% assembly of the time-stepping matrix
    K = spalloc(2*N_Omega, 2*N_Omega, 2*N_Omega*N_Omega);
    
        % first equation
        K(1:N_Omega,1:N_Omega) = M_bulk;
        K(1:N_Omega,(N_Omega+1):(2 * N_Omega)) = tau * sigma * S_bulk;
        
        % second equation
        K((N_Omega + 1):(2 * N_Omega),1:N_Omega) = epsilon * S_bulk;
        K((N_Omega + 1):(2 * N_Omega),(N_Omega + 1):(2 * N_Omega)) = - M_bulk;
        
    ML_bulk = assembly_ML_HN1(Nodes_bulk,Elements_bulk,N_Omega);
    
    % Anfangswerte
    sol_alpha = alpha_0;
    sol_beta = zeros(N_Omega,1);
    
    sol_n = ones(2*N_Omega,1);
    
    % get u_1 by using the method of order 1
    
    alpha_n = sol_alpha(:,1);
    
    v_1 = - M_bulk * alpha_n;
    v_2 = (1/epsilon) * ML_bulk * W_minus_bulk(alpha_n);
        
    b = [v_1; v_2]; % rhs
        
    % true solution for CH HN1
    %{
    fun = @(x) (K * x) + b + (1/epsilon) * f_nonlin_HN(ML_bulk, x, N_Omega) ; % final funtion to solve
    %}
    
    % Inhomogenous rhs for convergence plots test OG
    %
    i=1;
    %rhs1 = @(x) -exp(-tau) * (norm(x)^2 -1)^2 - 8 * sigma * exp(-tau) * (2*norm(x)^2 -1);
    %rhs1 = @(x) -exp(-tau) * (norm(x)^2 -1)^2 - 8 * sigma * exp(-tau) * (2*norm(x)^2 -1);
    rhs1 = @(x) exp(i*tau) * cos(pi*norm(x)^2) + sigma * 4*exp(i*tau) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2));
    v1 = zeros(N_Omega,1);
        
    for j = 1:N_Omega
        v1(j) = rhs1(Nodes_bulk(j,:));
    end
        
    v1 = tau * M_bulk * v1; % Mass lumping for approximation
        
    %rhs2 = @(x) - 8 * epsilon * exp(-tau) * (2*norm(x)^2 -1) + 1/epsilon * (exp(-3*tau) * (norm(x)^2-1)^6 - exp(-tau) * (norm(x)^2-1)^2) -exp(-tau) * (norm(x)^2 -1)^2;
    %rhs2 = @(x) - 8 * epsilon * exp(-tau) * (2*norm(x)^2 -1) + 1/epsilon * ((exp(-tau) * (norm(x)^2-1)^2 +1)^3 - exp(-tau) * (norm(x)^2-1)^2 -1) -exp(-tau) * (norm(x)^2 -1)^2 -1;
    rhs2 = @(x) epsilon * 4*exp(i*tau) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2)) - exp(i*tau) * cos(pi*norm(x)^2); %+ 1/epsilon * ((exp(i*tau) * cos(pi*norm(x)^2))^3 - exp(i*tau) * cos(pi*norm(x)^2));
    v2 = zeros(N_Omega,1);
    
    for j = 1:N_Omega
        v2(j) = rhs2(Nodes_bulk(j,:));
    end
        
    v2 = M_bulk * v2; % Mass lumping for approximation
    v = [v1;v2];
    fun = @(x) (K * x) + b + (1/epsilon) * f_nonlin_HN(ML_bulk, x, N_Omega) - v ; % final funtion to solve
    %}
        
    
    tol = 0.00000001;
    sol_n = newton_solver_HN1(ML_bulk, K, fun, N_Omega, epsilon, sol_n, tol);
        
    % update the solution vectors
    %{
    sol_alpha = [sol_alpha, sol_n(1:N_Omega)];
    sol_beta = [sol_beta, sol_n((N_Omega + 1):(2 * N_Omega))];
    %}
    
    % manipulate the solution vectors
    %
    temp = zeros(N_Omega,1);
    fctn = @(x,t) exp(t) * cos(pi*norm(x)^2);
    for i = 1:N_Omega
        temp(i) = fctn(Nodes_bulk(i,:),tau);
    end
    sol_alpha = [sol_alpha,temp];
    sol_beta = [sol_beta,temp];
    %}
    
    %% now order 2
    
    weights = getweights(0,1,nr);
    
    %% assembly of the time-stepping matrix
    K = spalloc(2*N_Omega, 2*N_Omega, 2*N_Omega*N_Omega);
    
        % first equation
        K(1:N_Omega,1:N_Omega) = M_bulk;
        K(1:N_Omega,(N_Omega+1):(2 * N_Omega)) = tau * sigma * (1/2) * S_bulk;
        
        % second equation
        K((N_Omega + 1):(2 * N_Omega),1:N_Omega) = epsilon * (1/2) * S_bulk;
        K((N_Omega + 1):(2 * N_Omega),(N_Omega + 1):(2 * N_Omega)) = - (1/2) * M_bulk;
        
    %% loop
    
    for i = 2:N
        
        i % Step
        
        alpha_n = sol_alpha(:,i);
        alpha_pn = sol_alpha(:,i-1);
        
        beta_n = sol_beta(:,i);
        
        v_1 = - M_bulk * alpha_n + tau * sigma * (1/2) * S_bulk * beta_n;
        v_2 = epsilon * (1/2) * S_bulk * alpha_n + (1/epsilon) * ((3/2) * M_bulk * W_minus_bulk(alpha_n) - (1/2) * M_bulk * W_minus_bulk(alpha_pn)) - (1/2) * M_bulk * beta_n;
        
        b = [v_1; v_2]; % rhs
        
        % true solution for CH HN2
        %{
        fun = @(x) (K * x) + b + f_nonlin_HN2(ML_bulk, x, N_Omega, epsilon, weights, alpha_n) ; % final funtion to solve
        %}
        
        % Inhomogenous rhs for convergence plots test OG
        %
        %rhs1 = @(x,t) exp(t) * (norm(x)^2 -1)^2 - 8 * sigma * exp(t) * (2*norm(x)^2 -1);
        %rhs1 = @(x) -exp(-i*tau) * (norm(x)^2 -1)^2 - 8 * sigma * exp(-i*tau) * (2*norm(x)^2 -1);
        rhs1 = @(x,t) exp(t) * cos(pi*norm(x)^2) + sigma * 4*pi*exp(t) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2));
        v1 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v1(j) = rhs1(Nodes_bulk(j,:),(i+0.5)*tau);
        end
        
        v1 = tau * M_bulk * v1; % Mass lumping for approximation
        
        %rhs2 = @(x,t) - 8 * epsilon * exp(t) * (2*norm(x)^2 -1) + 1/epsilon * (exp(3*t) * (norm(x)^2-1)^6 - exp(t) * (norm(x)^2-1)^2) -exp(t) * (norm(x)^2 -1)^2;
        %rhs2 = @(x) - 8 * epsilon * exp(-i*tau) * (2*norm(x)^2 -1) + 1/epsilon * ((exp(-i*tau) * (norm(x)^2-1)^2 +1)^3 - exp(-i*tau) * (norm(x)^2-1)^2 -1) -exp(-i*tau) * (norm(x)^2 -1)^2 -1;
        rhs2 = @(x,t) epsilon * 4*pi*exp(t) * (sin(pi*norm(x)^2) + pi*norm(x)^2 * cos(pi* norm(x)^2)) - exp(t) * cos(pi*norm(x)^2) + 1/epsilon * ((exp(t) * cos(pi*norm(x)^2))^3 - exp(t) * cos(pi*norm(x)^2));
        v2 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v2(j) = rhs2(Nodes_bulk(j,:),(i+0.5)*tau);
        end
        
        v2 = M_bulk * v2; % Mass lumping for approximation
        v = [v1;v2];
        fun = @(x) (K * x) + b - v + f_nonlin_HN2(M_bulk, x, N_Omega, epsilon, weights, alpha_n); % final funtion to solve
        %}
        
        % sol_n = fsolve(fun, sol_n);
        
        %
        tol = 0.000000001;
        sol_n = newton_solver_HN2(M_bulk, K, fun, N_Omega, epsilon, sol_n, tol, weights, alpha_n);
        %}
        
        %{
        sol_n = K\(-b+v);
        %}
        
        % update solution vector
        % alpha
        % sol((i * N_Omega + 1):((i + 1) * N_Omega),1) = sol_n(1:N_Omega);
         sol_alpha = [sol_alpha, sol_n(1:N_Omega)];
        % beta
        % sol(((N + i) * N_Omega + 1):((N + i + 1) * N_Omega),1) = sol_n((N_Omega + 1):(2 * N_Omega));
         sol_beta = [sol_beta, sol_n((N_Omega + 1):(2 * N_Omega))];
        
    end
end