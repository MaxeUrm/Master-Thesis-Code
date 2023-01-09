% function to assemble the mass lumping matrix

function [ML_bulk,ML_surf] = assembly_ML(Nodes_bulk,Elements_bulk,Nodes_surf,Elements_surf,N_Omega,N_Gamma)

    % Mass lumping Matrix in the bulk
    
    No_Elements = length(Elements_bulk);
    v = zeros(N_Omega,1);
    
    for i = 1:N_Omega
        for j = 1:No_Elements
            Nodes = Elements_bulk(j,:);
            if ismember(i,Nodes)  % check if Node i is in the element j
                x1 = Nodes_bulk(Elements_bulk(j,1),:);
                x2 = Nodes_bulk(Elements_bulk(j,2),:);
                x3 = Nodes_bulk(Elements_bulk(j,3),:);
        
                M = [x1;x2;x3];
        
                v(i) = v(i) + polyarea(M(:,1),M(:,2));
            end
        end
    end
    
    v = 1/3 * v;
    
    ML_bulk = diag(v);
    
    % Mass lumping Matrix on the surface
    
    No_Elements = length(Elements_surf);
    v = zeros(N_Gamma,1);
    
    for i = 1:N_Gamma
        for j = 1:No_Elements
            Nodes = Elements_surf(j,:);
            if ismember(i,Nodes)  % check if Node i is in the element j
                x1 = Nodes_surf(Elements_surf(j,1),:);
                x2 = Nodes_surf(Elements_surf(j,2),:);
                
                v(i) = v(i) + norm(x1-x2);
            end
        end
    end
    
    v = 1/2 * v;
    
    ML_surf = diag(v);
    
end