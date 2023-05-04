% Newton Method solver for HN2

function sol = newton_solver_HN2(ML_bulk, K, fun, N_Omega, epsilon, sol_n, tol, weights, alpha_n)
    
    sol = sol_n;
    
    n = length(weights);
    
    while norm(fun(sol)) > tol
        
        f_deriv = spalloc(2*N_Omega, 2*N_Omega, 2*N_Omega*N_Omega);
    
        alpha_next = sol(1:N_Omega);
        
        v1 = zeros(N_Omega,1);
        
        for j=1:n
            v1 = v1 + weights(j) * ((j-1)/(n-1)) * ML_bulk * W_plus_bulk_2((1-(j-1)/(n-1))*alpha_n + ((j-1)/(n-1))*alpha_next);
        end
        
        M1 = diag(v1);
        
        M1 = (1/epsilon) * M1;
        
        f_deriv((N_Omega +1):(2*N_Omega),1:N_Omega) = M1;
        
        f_deriv = f_deriv + K; % Jacobi Matrix of the function fun
        
        del_x = f_deriv\(-fun(sol));
        
        sol = sol + del_x;
        
        clear f_deriv;
    end

end