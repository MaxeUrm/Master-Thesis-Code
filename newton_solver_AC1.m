% Newton Method solver for AC

function sol = newton_solver_AC1(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, tau, delta, sol_n, tol)
    
    sol = sol_n;
    
    while norm(fun(sol)) > tol
        
        f_deriv = spalloc(2*(N_Omega + N_Gamma), 2*(N_Omega + N_Gamma), 8*(N_Omega + N_Gamma));
        
        x1 = sol(1:N_Omega);
        x2 = sol((2*N_Omega + N_Gamma + 1):end);
        
        M1 = diag(W_plus_bulk_2(x1));
        M2 = diag(W_plus_surf_2(x2));
        
        M1 = (1/epsilon) * ML_bulk * M1;
        M2 = (tau/delta) * ML_surf * M2;
        
        f_deriv((N_Omega +1):(2*N_Omega),1:N_Omega) = M1;
        f_deriv((2*N_Omega +1):(2*N_Omega + N_Gamma),(2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma)) = M2;
        
        f_deriv = f_deriv + K; % Jacobi Matrix of the function fun
        
        del_x = f_deriv\(-fun(sol));
        
        sol = sol + del_x;
        
        clear f_deriv;
    end

end