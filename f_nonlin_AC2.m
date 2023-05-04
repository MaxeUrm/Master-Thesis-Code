% Nonlinear function for AC2

function sol = f_nonlin_AC2(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, tau, delta, weights, alpha_n, xi_n)
    
    alpha_next = x(1:N_Omega);
    xi_next = x((2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma));
    
    n = length(weights);
    
    v1 = zeros(N_Omega,1);
    v2 = zeros(N_Gamma,1);
    
    for j = 1:n
        vec = zeros(N_Omega,1);
        
        vec = W_plus_bulk((1-((j-1)/(n-1)))*alpha_n + ((j-1)/(n-1))*alpha_next);
        
        v1 = v1 + weights(j) * vec;
        
        vec = zeros(N_Gamma,1);
        
        vec = W_plus_surf((1-((j-1)/(n-1)))*xi_n + ((j-1)/(n-1))*xi_next);
        
        v2 = v2 + weights(j) * vec;
    end
    
    v1 = 1/epsilon * ML_bulk * v1;
    v2 = tau * (1/delta) * ML_surf * v2;
    
    
    sol = [zeros(N_Omega,1); v1; v2; zeros(N_Gamma,1)];
    
end