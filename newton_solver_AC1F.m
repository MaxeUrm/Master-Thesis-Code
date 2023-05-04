% Newton Method solver for AC1F

function sol = newton_solver_AC1F(ML_bulk, ML_surf, K, fun, N_Omega, N_Gamma, epsilon, tau, delta, sol_n, tol, k, l)
    
    sol = sol_n;
    
    while norm(fun(sol)) > tol
        
        f_deriv = spalloc((2*k*N_Omega + (l+1)*N_Gamma), (2*k*N_Omega + (l+1)*N_Gamma), 4*(k+l)*(N_Omega + N_Gamma));
    
        for i = 1:k
            v1 = sol(((i-1)*N_Omega +1):(i*N_Omega),1); % alpha n+i/k
            v1 = W_plus_bulk_2(v1);
            M1 = diag(v1);
            M1 = (1/epsilon) * ML_bulk .* M1;
            f_deriv(((i-1+k)*N_Omega+1):((i+k)*N_Omega),((i-1)*N_Omega +1):(i*N_Omega)) = M1;
            if i > 1
                v0 = sol(((i-2)*N_Omega +1):((i-1)*N_Omega)); % alpha n+(i-1)/k
                v0 = W_minus_bulk_2(v0);
                M0 = diag(v0);
                M0 = (1/epsilon) * ML_bulk .* M0;
                f_deriv(((i-1+k)*N_Omega+1):((i+k)*N_Omega),((i-2)*N_Omega +1):((i-1)*N_Omega)) = M0;
            end
        end
        
        for i = 1:l
            v1 = sol((2*k*N_Omega + i*N_Gamma +1):(2*k*N_Omega + (i+1)*N_Gamma),1); % xi n+i/l
            v1 = W_plus_surf_2(v1);
            M1 = diag(v1);
            M1 = (tau/(delta*l)) * ML_surf .* M1;
            f_deriv((2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma),(2*k*N_Omega + i*N_Gamma +1):(2*k*N_Omega + (i+1)*N_Gamma)) = M1;
            if i > 1
                v0 = sol((2*k*N_Omega + (i-1)*N_Gamma +1):(2*k*N_Omega + i*N_Gamma),1); % xi n+(i-1)/l
                v0 = W_minus_surf_2(v0);
                M0 = diag(v0);
                M0 = (tau/(delta*l)) * ML_surf .* M0;
                f_deriv((2*k*N_Omega + (i-1)*N_Gamma + 1):(2*k*N_Omega + i*N_Gamma),(2*k*N_Omega + (i-1)*N_Gamma +1):(2*k*N_Omega + i*N_Gamma)) = M0;
            end
        end
        
        f_deriv = f_deriv + K; % Jacobi Matrix of the function fun
        
        del_x = f_deriv\(-fun(sol));
        
        sol = sol + del_x;
        
        clear f_deriv;
    end

end