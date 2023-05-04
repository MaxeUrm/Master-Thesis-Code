%% Time-stepping scheme slover for CH-System with Allen-Cahn-Type Boundary

% Order 1

function [sol_alpha,sol_beta,sol_gamma,sol_xi] = time_stepping_AC1FF(S_bulk,M_bulk,Nodes_bulk,Elements_bulk,S_surf,M_surf,B,B_I, Nodes_surf, Elements_surf, alpha_0, xi_0, T, N, k, l, epsilon, sigma, delta, kappa);

    %% dof in the bulk
    N_Omega = length(Nodes_bulk);
    % dof on the surface
    N_Gamma = length(Nodes_surf);
    
    %% tau
    tau = T/N;
    
    %% assembly of the time-stepping matrix
    K = spalloc((2*k*N_Omega + (l+1)*N_Gamma), (2*k*N_Omega + (l+1)*N_Gamma), 4*(k+l)*(N_Omega + N_Gamma));
    
        % first equation
        for i = 1:k
            K(((i-1)*N_Omega +1):(i * N_Omega),((i-1)*N_Omega +1):(i * N_Omega)) = M_bulk;
            K(((i-1)*N_Omega +1):(i * N_Omega),((i-1+k)*N_Omega +1):((i+k) * N_Omega)) = (tau/k) * sigma * S_bulk;
            if i > 1
                K(((i-1)*N_Omega +1):(i * N_Omega),((i-2)*N_Omega +1):((i-1) * N_Omega)) = -M_bulk;
            end
        end
        
        % second equation
        for i = 1:k
            K(((i-1+k)*N_Omega +1):((i+k) * N_Omega),((i-1)*N_Omega +1):(i * N_Omega)) = epsilon * S_bulk;
            K(((i-1+k)*N_Omega +1):((i+k) * N_Omega),((i-1+k)*N_Omega +1):((i+k) * N_Omega)) = -M_bulk;
            K(((i-1+k)*N_Omega +1):((i+k) * N_Omega),(2*k*N_Omega +1):(2*k*N_Omega + N_Gamma)) = -epsilon * B;
        end
        
        % thrid equation
        for i = 1:l
            K((2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma),(2*k*N_Omega + i*N_Gamma + 1):(2*k*N_Omega + (i+1)*N_Gamma)) = M_surf + (tau/l) * delta * kappa * S_surf;
            K((2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma),(2*k*N_Omega +1):(2*k*N_Omega + N_Gamma)) = epsilon * (tau/l) * M_surf;
            if i > 1
                K((2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma),(2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma)) = - M_surf;
            end
        end
        
        % fourth equation
        K((2*k*N_Omega + l*N_Gamma +1):(2*k*N_Omega + (l+1)*N_Gamma),((k-1)*N_Omega +1):(k * N_Omega)) = B_I;
        K((2*k*N_Omega + l*N_Gamma +1):(2*k*N_Omega + (l+1)*N_Gamma),(2*k*N_Omega + l*N_Gamma + 1):(2*k*N_Omega + (l+1)*N_Gamma)) = - eye(N_Gamma);
        
    [ML_bulk,ML_surf] = assembly_ML(Nodes_bulk,Elements_bulk,Nodes_surf,Elements_surf,N_Omega,N_Gamma);
    
    
    %% loop
    
    % Anfangswerte
    sol_alpha = alpha_0;
    sol_beta = zeros(N_Omega,1);
    sol_gamma = zeros(N_Gamma,1);
    sol_xi = xi_0;
    
    sol_n = ones(2*k*N_Omega + (l+1)*N_Gamma,1);
    
    for i = 1:N
        
        i % Step
        
        % construct rhs out of prevous step
        
        alpha_n = sol_alpha(:,(i-1)*k +1);
        xi_n = sol_xi(:,(i-1)*l +1);
        
        % these are needed for the prevous timestep
        v_1 = - M_bulk * alpha_n;
        v_2 = (1/epsilon) * ML_bulk * W_minus_bulk(alpha_n);
        v_3 = - M_surf * xi_n + (tau/(delta*l)) * ML_surf * W_minus_surf(xi_n);
        v_4 = zeros(l*N_Gamma,1);
        
        b = [v_1; zeros((k-1)*N_Omega,1); v_2; zeros((k-1)*N_Omega,1); v_3; v_4]; % rhs
        
        % true solution to CH AC1FF
        %{
        fun = @(x) K * x + b + f_nonlin_AC1F(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta, k, l); % final funtion to solve
        %}
        
        % error improvement section
        %
        rhs1 = @(x,t) -exp(-t) * (norm(x)^2 + x(1)*x(2)) - 8 * sigma * exp(-t) * (2*norm(x)^2 -1);
        
        % rhs1 = @(x,t) -exp(-t);
        
        v1 = zeros(k*N_Omega,1);
        for r = 1:k
            temp = zeros(N_Omega,1);
            for j = 1:N_Omega
                temp(j) = rhs1(Nodes_bulk(j,:),(i+(r/k)-1)*tau);
            end
        
            temp = (tau/k) * ML_bulk * temp; % Mass lumping for approximation
            
            v1(((r-1)*N_Omega +1):r*N_Omega) = temp;
        end 
        
        rhs2 = @(x,t) - 4 * epsilon * exp(-t) + 1/epsilon * (exp(-3*t) * (norm(x)^2 + x(1)*x(2))^3 - exp(-t) * (norm(x)^2 + x(1)*x(2))) - exp(-t) * (norm(x)^2 -1)^2;
        
        % rhs2 = @(x,t) 1/epsilon * (exp(-3*t) - exp(-t));
        
        v2 = zeros(k*N_Omega,1);
        for r = 1:k
            temp = zeros(N_Omega,1);
            for j = 1:N_Omega
                temp(j) = rhs2(Nodes_bulk(j,:),(i+(r/k)-1)*tau);
            end
        
            temp = ML_bulk * temp; % Mass lumping for approximation
            
            v2(((r-1)*N_Omega +1):r*N_Omega) = temp;
        end 
        
        rhs3 = @(x,t) -exp(-t) * (1+x(1)*x(2)) + 4 * delta * kappa * exp(-t) * x(1)*x(2) + 1/delta * (exp(-3*t) * (1+x(1)*x(2))^3 - exp(-t) * (1+x(1)*x(2))) + 2 * epsilon * exp(-t) * (1+x(1)*x(2));
        
        % rhs3 = @(x,t) -exp(-t) + 1/delta * (exp(-3*t) - exp(-t));
        
        v3 = zeros(l*N_Gamma,1);
        for r = 1:l
            temp = zeros(N_Gamma,1);
            for j = 1:N_Gamma
                temp(j) = rhs3(Nodes_surf(j,:),(i+(r/l)-1)*tau);
            end
        
            temp = (tau/l) * ML_surf * temp; % Mass lumping for approximation
            
            v3(((r-1)*N_Gamma +1):r*N_Gamma) = temp;
        end 
        
        
        v = [v1;v2;v3;zeros(N_Gamma,1)];
        fun = @(x) (K * x) + b + f_nonlin_AC1F(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta, k, l) - v ; % final funtion to solve
        %}
        
        % sol_n = fsolve(fun, sol_n); % solving...
        
        
        tol = 0.000000001;
        sol_n = newton_solver_AC1F(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, tau, delta, sol_n, tol, k, l);
        
        
        % update solution vector
        % alpha
        for j = 1:k
            sol_alpha = [sol_alpha,sol_n(((j-1)*N_Omega +1):(j*N_Omega),1)];
        end
        % beta
        for j = 1:k
            sol_beta = [sol_beta,sol_n(((j-1+k)*N_Omega +1):((j+k)*N_Omega),1)];
        end
        % gamma
        sol_gamma = [sol_gamma,sol_n((2*k*N_Omega +1):(2*k*N_Omega + N_Gamma),1)];
        % xi
        for j = 1:l
            sol_xi = [sol_xi,sol_n((2*k*N_Omega + j*N_Gamma +1):(2*k*N_Omega + (j+1)*N_Gamma),1)];
        end
        
    end
end