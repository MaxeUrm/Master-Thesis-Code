% Nonlinear function for GMS2

function sol = f_nonlin_GMS2(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta, weights, alpha_n, xi_n)
    
    alpha_next = x(1:N_Omega);
    xi_next = x((2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma));
    
    n = length(weights);
    
    v1 = zeros(N_Omega,1);
    v2 = zeros(N_Gamma,1);
    
    for j = 1:n
        vec = zeros(N_Omega,1);
        
        for i = 1:N_Omega
            vec(i) = W_plus_bulk((1-((j-1)/(n-1)))*alpha_n(i) + ((j-1)/(n-1))*alpha_next(i));
        end
        
        v1 = v1 + weights(j) * ML_bulk * vec;
        
        vec = zeros(N_Gamma,1);
        
        for i = 1:N_Gamma
            vec(i) = W_plus_surf((1-((j-1)/(n-1)))*xi_n(i) + ((j-1)/(n-1))*xi_next(i));
        end
        
        v2 = v2 + weights(j) * ML_surf * vec;
    end
    
    v1 = 1/epsilon * v1;
    v2 = 1/delta * v2;
    
    
    sol = [zeros(N_Omega,1); v1; zeros(N_Gamma,1); v2; zeros(N_Gamma,1); zeros(N_Gamma,1)];
    
end