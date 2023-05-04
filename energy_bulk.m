% Energy functional bulk

function sol = energy_bulk(val, epsilon, Nodes_bulk, Elements_bulk, S_bulk)
    % first part
    N_Omega = length(val);
    % M = val * transpose(val);
    % M = M .* S_bulk;
    
    % sol = (epsilon/2) * ones(1,N_Omega) * M * ones(N_Omega,1); % add all values
    sol = (epsilon/2) * transpose(val) * S_bulk * val;
    
    % second part
    v = W(val);
    
    sol = sol + (1/epsilon) * integral_approx2d(v,Nodes_bulk,Elements_bulk);
    
end