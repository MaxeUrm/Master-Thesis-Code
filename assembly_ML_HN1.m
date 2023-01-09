% function to assemble the mass lumping matrix for HN1

function ML_bulk = assembly_ML_HN1(Nodes_bulk,Elements_bulk,N_Omega)

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
        
                v(i) = polyarea(M(:,1),M(:,2));
            end
        end
    end
    
    v = 1/3 * v;
    
    ML_bulk = diag(v);
    
end