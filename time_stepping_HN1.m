%% Time-stepping scheme slover for CH-System with homogenous Neumann boundary values
% values

% Order 1

function [sol_alpha,sol_beta] = time_stepping_HN1(S_bulk,M_bulk,Nodes_bulk, Elements_bulk, alpha_0, T, N, epsilon, sigma)
    
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
    
    
    %% loop
    
    % Anfangswerte
    % sol(1:N_Omega,1) = alpha_0;
    sol_n = ones(2*N_Omega,1);
    
    sol_alpha = alpha_0;
    sol_beta = zeros(N_Omega,1);
    
    for i = 1:N
        
        %i % Step
        
        % construct rhs out of prevous step
        
        % alpha_n = sol(((i-1) * N_Omega + 1):(i * N_Omega),1);
        alpha_n = sol_alpha(:,i);
        
        v_1 = - M_bulk * alpha_n;
        v_2 = (1/epsilon) * ML_bulk * W_minus_bulk(alpha_n);
        
        b = [v_1; v_2]; % rhs
        
        % true solution for CH HN1
        %{
        fun = @(x) (K * x) + b + (1/epsilon) * f_nonlin_HN(ML_bulk, x, N_Omega) ; % final funtion to solve
        %}
        
        % Inhomogenous rhs for convergence plots test OG
        %
        rhs1 = @(x) -exp(-i*tau) * (norm(x)^2 -1)^2 - 8 * sigma * exp(-i*tau) * (2*norm(x)^2 -1);
        v1 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v1(j) = rhs1(Nodes_bulk(j,:));
        end
        
        v1 = tau * ML_bulk * v1; % Mass lumping for approximation
        
        rhs2 = @(x) - 8 * epsilon * exp(-i*tau) * (2*norm(x)^2 -1) + 1/epsilon * (exp(-3*i*tau) * (norm(x)^2-1)^6 - exp(-i*tau) * (norm(x)^2-1)^2) -exp(-i*tau) * (norm(x)^2 -1)^2;
        v2 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v2(j) = rhs2(Nodes_bulk(j,:));
        end
        
        v2 = ML_bulk * v2; % Mass lumping for approximation
        v = [v1;v2];
        fun = @(x) (K * x) + b + (1/epsilon) * f_nonlin_HN(ML_bulk, x, N_Omega) - v ; % final funtion to solve
        %}
        
        % sol_n = fsolve(fun, sol_n);
        
        tol = 0.000001;
        sol_n = newton_solver_HN1(ML_bulk, K, fun, N_Omega, epsilon, sol_n, tol);
        
        % update solution vector
        % alpha
        % sol((i * N_Omega + 1):((i + 1) * N_Omega),1) = sol_n(1:N_Omega);
         sol_alpha = [sol_alpha, sol_n(1:N_Omega)];
        % beta
        % sol(((N + i) * N_Omega + 1):((N + i + 1) * N_Omega),1) = sol_n((N_Omega + 1):(2 * N_Omega));
         sol_beta = [sol_beta, sol_n((N_Omega + 1):(2 * N_Omega))];
        
    end
end