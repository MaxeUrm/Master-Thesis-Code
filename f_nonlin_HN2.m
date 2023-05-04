% Nonlinear function

function sol = f_nonlin_HN2(ML_bulk, x, N_Omega, epsilon, weights, alpha_n)

    alpha_next = x(1:N_Omega);
    
    n = length(weights);
    
    v1 = zeros(N_Omega,1);
    
    for j = 1:n
        vec = zeros(N_Omega,1);
        
        for i = 1:N_Omega
            vec(i) = W_plus_bulk((1-((j-1)/(n-1)))*alpha_n(i) + ((j-1)/(n-1))*alpha_next(i));
        end
        
        v1 = v1 + weights(j) * ML_bulk * vec;
        
    end
    
    v1 = 1/epsilon * v1;
    
    sol = [zeros(N_Omega,1); v1];
    
end