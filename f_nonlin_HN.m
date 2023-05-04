% Nonlinear function

function sol = f_nonlin_HN(ML_bulk, x, N_Omega)
    v1 = x(1:N_Omega);
    
    v1 = W_plus_bulk(v1);
    
    v1 = ML_bulk * v1;
    
    sol = [zeros(N_Omega,1); v1];
    
end