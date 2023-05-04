%% Time-stepping scheme slover for CH-System with Liu and Wu Boundary values

% Order 1

function [sol_alpha,sol_beta,sol_gamma,sol_xi,sol_zeta,sol_eta] = time_stepping_GMS1(S_bulk,M_bulk,S_surf,M_surf,Nodes_bulk, Elements_bulk, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, epsilon, sigma, delta, kappa)
    
    %% dof in the bulk
    N_Omega = length(Nodes_bulk);
    % dof on the surface
    N_Gamma = length(Nodes_surf);
    
    %% tau
    tau = T/N;
    
    %% assembly of the time-stepping matrix
    K = spalloc((2*N_Omega + 4*N_Gamma), (2*N_Omega + 4*N_Gamma), 8*(N_Omega + N_Gamma));
    
    B = assembly_B(Nodes_bulk, Nodes_surf);
    
        % first equation
        K(1:N_Omega,1:N_Omega) = M_bulk;
        K(1:N_Omega,(N_Omega+1):(2 * N_Omega)) = tau * sigma * S_bulk;
        K(1:N_Omega,(2 * N_Omega + 3 * N_Gamma + 1):(2 * N_Omega + 4 * N_Gamma)) = -tau * transpose(B) * M_surf;
        
        % second equation
        K((N_Omega + 1):(2 * N_Omega),1:N_Omega) = epsilon * S_bulk;
        K((N_Omega + 1):(2 * N_Omega),(N_Omega + 1):(2 * N_Omega)) = -M_bulk;
        K((N_Omega + 1):(2 * N_Omega),(2 * N_Omega + 1):(2 * N_Omega + N_Gamma)) = -epsilon * transpose(B) * M_surf;
        
        % third equation
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = M_surf;
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma)) = tau * S_surf;
        K((2 * N_Omega + 1):(2 * N_Omega + N_Gamma),(2 * N_Omega + 3 * N_Gamma + 1):(2 * N_Omega + 4 * N_Gamma)) = tau * M_surf;
        
        % fourth equation
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + 1):(2 * N_Omega + N_Gamma)) = epsilon * M_surf;
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = delta * kappa * S_surf;
        K((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma)) = - M_surf;
        
        % fifth equation
        K((2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma),1:N_Omega) = B;
        K((2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = - diag(ones(N_Gamma,1));
        
        % sixth equation
        K((2 * N_Omega + 3 * N_Gamma + 1):(2 * N_Omega + 4 * N_Gamma),(N_Omega + 1):(2 * N_Omega)) = B;
        K((2 * N_Omega + 3 * N_Gamma + 1):(2 * N_Omega + 4 * N_Gamma),(2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma)) = - diag(ones(N_Gamma,1));
        
    [ML_bulk,ML_surf] = assembly_ML(Nodes_bulk,Elements_bulk,Nodes_surf,Elements_surf,N_Omega,N_Gamma);
    
    %% loop
    
    % Anfangswerte
    sol_alpha = alpha_0;
    sol_beta = zeros(N_Omega,1);
    sol_gamma = zeros(N_Gamma,1);
    sol_xi = xi_0;
    sol_zeta = zeros(N_Gamma,1);
    sol_eta = zeros(N_Gamma,1);
    
    sol_n = ones(2*N_Omega+4*N_Gamma,1);
    
    for i = 1:N
        
        i % Step
        
        % construct rhs out of prevous step
        
        alpha_n = sol_alpha(:,i);
        xi_n = sol_xi(:,i);
        
        v_1 = - M_bulk * alpha_n;
        v_2 = (1/epsilon) * ML_bulk * W_minus_bulk(alpha_n);
        v_3 = - M_surf * xi_n;
        v_4 = (1/delta) * ML_surf * W_minus_surf(xi_n);
        v_5 = zeros(N_Gamma,1);
        
        b = [v_1; v_2; v_3; v_4; v_5; v_5]; % rhs
        
        % true solution forCH GMS1
        %{
        fun = @(x) K * x + b + f_nonlin_GMS1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta); % final funtion to solve
        %}
        
        % Inhomogenous rhs for convergence plot tests OG
        %
        rhs1 = @(x) -exp(-i*tau) * norm(x)^2 - 4 * sigma * exp(-i*tau);
        v1 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v1(j) = rhs1(Nodes_bulk(j,:));
        end
        
        v1 = tau * ML_bulk * v1; % Mass lumping for approximation
        
        rhs2 = @(x) - 4 * epsilon * exp(-i*tau) + 1/epsilon * (exp(-3*i*tau) * norm(x)^6 - exp(-i*tau) * norm(x)^2) -exp(-i*tau) * norm(x)^2;
        v2 = zeros(N_Omega,1);
        
        for j = 1:N_Omega
            v2(j) = rhs2(Nodes_bulk(j,:));
        end
        
        v2 = ML_bulk * v2; % Mass lumping for approximation
        
        rhs3 = @(x) - exp(-i*tau) + 2 * sigma * exp(-i*tau);
        v3 = zeros(N_Gamma,1);
        
        for j = 1:N_Gamma
            v3(j) = rhs3(Nodes_surf(j,:));
        end
        v3 = tau * ML_surf * v3;
        
        
        rhs4 = @(x) - 4 * delta * kappa * exp(-i*tau) * 0 + 1/delta * (exp(-3*i*tau) - exp(-i*tau)) + 2 * epsilon * exp(-i*tau) - exp(-i*tau);
        v4 = zeros(N_Gamma,1);
        
        for j = 1:N_Gamma
            v4(j) = rhs4(Nodes_surf(j,:));
        end
        v4 = ML_surf * v4;
        
        v = [v1;v2;v3;v4;zeros(2*N_Gamma,1)];
        fun = @(x) (K * x) + b + f_nonlin_GMS1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta) - v ; % final funtion to solve
        %}
        
        % sol_n = fsolve(fun, sol_n);
        
        tol = 0.000001;
        sol_n = newton_solver_GMS1(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, delta, sol_n, tol);
        
        % update solution vector
        % alpha
        sol_alpha = [sol_alpha,sol_n(1:N_Omega)];
        % beta
        sol_beta = [sol_beta,sol_n((N_Omega + 1):(2 * N_Omega))];
        % gamma
        sol_gamma = [sol_gamma,sol_n((2 * N_Omega + 1):(2 * N_Omega + N_Gamma))];
        % xi
        sol_xi = [sol_xi,sol_n((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma))];
        % zeta
        sol_zeta = [sol_zeta,sol_n((2 * N_Omega + 2 * N_Gamma + 1):(2 * N_Omega + 3 * N_Gamma))];
        % eta
        sol_eta = [sol_eta,sol_n((2 * N_Omega + 3 * N_Gamma + 1):(2 * N_Omega + 4 * N_Gamma))];
        
        
    end
end