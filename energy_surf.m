% Energy functional surface

function sol = energy_surf(val, delta, kappa, Nodes_surf, Elements_surf, S_surf)
    % first part
    N_Gamma = length(val);
    M = val * transpose(val);
    M = M .* S_surf;
    
    sol = (delta*kappa/2) * ones(1,N_Gamma) * M * ones(N_Gamma,1); % add all components
    
    % second part
    v = W_Gamma(val);
    
    sol = sol + (1/delta) * integral_approx1d(v,Nodes_surf,Elements_surf);
    
end