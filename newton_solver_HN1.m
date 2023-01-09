% Newton Method solver for HN

function sol = newton_solver_HN1(ML_bulk, K, fun, N_Omega, epsilon, sol_n, tol)
    
    sol = sol_n;
    
    while norm(fun(sol)) > tol
        
        f_deriv = spalloc(2*N_Omega, 2*N_Omega, 2*N_Omega*N_Omega);
    
        x1 = sol(1:N_Omega);
        
        M1 = diag(W_plus_bulk_2(x1));
        
        M1 = (1/epsilon) * ML_bulk * M1;
        
        f_deriv((N_Omega +1):(2*N_Omega),1:N_Omega) = M1;
        f_deriv = f_deriv + K; % Jacobi Matrix of the function fun
        
        del_x = f_deriv\(-fun(sol));
        
        sol = sol + del_x;
        
        clear f_deriv;
    end

end