% this funcion gives an interpolation matrix to fullfill the connection
% term and to get initial values for the boundary

function B_I = assembly_interpolmatrix(Bdry, Nodes_surf, N_Omega, N_bdry, N_Gamma)

    B_I = spalloc(N_Gamma, N_Omega,2*N_bdry);
    N = length(Bdry(:,1));
    
    for i = 1:N_Gamma
        k = 1;
        while ~(i == Bdry(k,4))
            k = k+1; % get right node
        end
        
        if Bdry(k,5) == 1 % nodes are identical
            B_I(i,Bdry(k,3)) = 1;
        else
            
            n = k+1; % get next node with bulk origin
            if k+1 == N+1
            	n = 1;
            end
            
            while Bdry(n,5) == 0
                n = n+1;
                if n == N+1
                    n = 1;
                end
            end
            
            p = k-1; % get previous node with bulk origin
            if k-1 == 0
            	p = N;
            end
            
            while Bdry(p,5) == 0
                p = p-1;
                if p == 0
                    p = N;
                end
            end
            
            B_I(i,Bdry(p,3)) = norm(Bdry(p,1:2)-Bdry(k,1:2))/(norm(Bdry(p,1:2)-Bdry(k,1:2))+norm(Bdry(n,1:2)-Bdry(k,1:2)));
            
            B_I(i,Bdry(n,3)) = norm(Bdry(n,1:2)-Bdry(k,1:2))/(norm(Bdry(p,1:2)-Bdry(k,1:2))+norm(Bdry(n,1:2)-Bdry(k,1:2)));
        end
    end

end