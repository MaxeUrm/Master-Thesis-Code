% Assembly of the matrix that results from the trace operator

% This requires that surface nodes are also bulk nodes

function B = assembly_B (Nodes_bulk, Nodes_surf)
    
    % dof in the bulk
    N_Omega = length(Nodes_bulk);
    % dof on the surface
    N_Gamma = length(Nodes_surf);
    
    % introducing B
    B = spalloc(N_Gamma, N_Omega, N_Gamma); % Dimensions are N_Gamma x N_Omega with max N_Gamma ones
    
    for i = 1:N_Omega
        for j = 1:N_Gamma
            if Nodes_bulk(i,:) == Nodes_surf(j,:)
                B(j,i) = 1;
            end
        end
    end
    
end