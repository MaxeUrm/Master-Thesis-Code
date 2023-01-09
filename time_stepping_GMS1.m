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
        
    Z = spalloc((2*N_Omega + 4*N_Gamma), (2*N_Omega + 4*N_Gamma), 8*(N_Omega + N_Gamma));
    
    [ML_bulk,ML_surf] = assembly_ML(Nodes_bulk,Elements_bulk,Nodes_surf,Elements_surf,N_Omega,N_Gamma);
    
        Z(1:N_Omega,1:N_Omega) = M_bulk;
        Z((2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma),(2 * N_Omega + N_Gamma + 1):(2 * N_Omega + 2 * N_Gamma)) = M_surf;
    norm(Z,"fro")/norm(K,"fro")
    
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
        
        fun = @(x) K * x + b + f_nonlin_GMS1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta); % final funtion to solve
        
        % sol_n = fsolve(fun, sol_n);
        
        tol = 0.000000001;
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