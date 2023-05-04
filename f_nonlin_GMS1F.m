% Nonlinear function for GMS1F

function sol = f_nonlin_GMS1F(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta, k, l)
    
    sol = zeros(2*k*N_Omega + 2*(l+1)*N_Gamma,1);
    % second equation inlinearity
    for i = 1:k
        v1 = x(((i-1)*N_Omega +1):(i*N_Omega),1); % alpha n+i/k        
        v1 = W_plus_bulk(v1);
        v1 = (1/epsilon) * ML_bulk * v1;
        sol(((i-1+k)*N_Omega +1):((i+k)*N_Omega)) = sol(((i-1+k)*N_Omega +1):((i+k)*N_Omega)) + v1;
        if i > 1
            v0 = x(((i-2)*N_Omega +1):((i-1)*N_Omega)); % alpha n+(i-1)/k
            v0 = W_minus_bulk(v0);
            v0 = (1/epsilon) * ML_bulk * v0;
            sol(((i-1+k)*N_Omega +1):((i+k)*N_Omega),1) = sol(((i-1+k)*N_Omega +1):((i+k)*N_Omega),1) + v0;
        end
    end
    
    % fourth equation inlinearity
    for i = 1:l
        v1 = x((2*k*N_Omega + i*N_Gamma +1):(2*k*N_Omega + (i+1)*N_Gamma),1); % xi n+i/l
        v1 = W_plus_surf(v1);
        v1 = (1/delta) * ML_surf * v1;
        sol((2*k*N_Omega + (l+i-1)*N_Gamma + 1):(2*k*N_Omega + (l+i)*N_Gamma),1) = sol((2*k*N_Omega + (l+i-1)*N_Gamma + 1):(2*k*N_Omega + (l+i)*N_Gamma),1) + v1;
        if i > 1
            v0 = x((2*k*N_Omega + (i-1)*N_Gamma +1):(2*k*N_Omega + i*N_Gamma),1); % xi n+(i-1)/l
            v0 = W_minus_surf(v0);
            v0 = (1/delta) * ML_surf * v0;
            sol((2*k*N_Omega + (l+i-1)*N_Gamma + 1):(2*k*N_Omega + (l+i)*N_Gamma),1) = sol((2*k*N_Omega + (l+i-1)*N_Gamma + 1):(2*k*N_Omega + (l+i)*N_Gamma),1) + v0;
        end
    end
    

end