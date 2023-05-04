% Newton Method solver for GMS

function sol = newton_solver_GMS2(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, delta, sol_n, tol, weights, alpha_n, xi_n)
    
    sol = sol_n;
    
    n = length(weights);
    
    while norm(fun(sol)) > tol
        
        f_deriv = spalloc((2*N_Omega + 4*N_Gamma), (2*N_Omega + 4*N_Gamma), 10*(N_Omega + N_Gamma));
    
        alpha_next = sol(1:N_Omega);
        xi_next = sol((2*N_Omega + N_Gamma + 1):(2*N_Omega + 2*N_Gamma));
        
        v1 = zeros(N_Omega,1);
        
        for j=1:n
            v1 = v1 + weights(j) * ((j-1)/(n-1)) * ML_bulk * W_plus_bulk_2((1-(j-1)/(n-1))*alpha_n + ((j-1)/(n-1))*alpha_next);
        end
        
        v2 = zeros(N_Gamma,1);
        
        for j=1:n
            v2 = v2 + weights(j) * ((j-1)/(n-1)) * ML_surf * W_plus_surf_2((1-(j-1)/(n-1))*xi_n + ((j-1)/(n-1))*xi_next);
        end
        
        M1 = diag(v1);
        M2 = diag(v2);
        
        M1 = (1/epsilon) * M1;
        M2 = (1/delta) * M2;
        
        f_deriv((N_Omega +1):(2*N_Omega),1:N_Omega) = M1;
        f_deriv((2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma),(2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma)) = M2;
        
        f_deriv = f_deriv + K; % Jacobi Matrix of the function fun
        
        del_x = f_deriv\(-fun(sol));
        
        sol = sol + del_x;
        
        clear f_deriv;
    end

end