% Nonlinear function for LW

function sol = f_nonlin_GMS1(ML_bulk, ML_surf, x, N_Omega, N_Gamma, epsilon, delta)
    v1 = x(1:N_Omega);
    v2 = x((2*N_Omega + N_Gamma +1):(2*N_Omega + 2*N_Gamma));
    
    v1 = W_plus_bulk(v1);
    v2 = W_plus_bulk(v2);
    
    v1 = (1/epsilon) * ML_bulk * v1;
    v2 = (1/delta) * ML_surf * v2;
    
    sol = [zeros(N_Omega,1); v1; zeros(N_Gamma,1); v2; zeros(N_Gamma,1); zeros(N_Gamma,1)];
    
end