% This funktion builds the weakly defined Matrix B for different
% discretizations on the surface and on the boundary

function [B,Bdry] = assembly_tracematrix(Nodes_bulk,Bdrynodes_bulk,Nodes_surf,N_Omega,N_bdry,N_Gamma,position_array)

    B = spalloc(N_Omega,N_Gamma,N_bdry); % final matrix
    
    N = N_Gamma+N_bdry;
    
    Bdry = zeros(N,7); % collect the boundary nodes in order to intersect intervalls later
    
    for i = 1:N_bdry
        Bdry(i,1:2) = Bdrynodes_bulk(i,:);
        Bdry(i,3) = position_array(i);
        Bdry(i,5) = 1; % 1 indicates Node with inner origin
        Bdry(i,6) = 0;
        Bdry(i,7) = atan2(Bdrynodes_bulk(i,2),Bdrynodes_bulk(i,1)); % get angle with x axis
    end
    
    for i = 1:N_Gamma
        Bdry(i+N_bdry,1:2) = Nodes_surf(i,:);
        Bdry(i+N_bdry,4) = i;
        Bdry(i+N_bdry,5) = 0;
        Bdry(i+N_bdry,6) = 1; % 1 indicates Node with Gamma origin
        Bdry(i+N_bdry,7) = atan2(Nodes_surf(i,2),Nodes_surf(i,1)); % get angle with x axis (to sort like this only works in certain domains)
    end
    
    Bdry = sortrows(Bdry, 7); % Nodes are sortet counterclockwise on the unified boundary
    
    temp = Bdry;
    nr = 0;
    
    for i = 1:(N-1) % check for double entries i.e. if Gamma node is a Omega node
        if Bdry(i,7) == Bdry(i+1,7)
            temp(i-nr,3) = Bdry(i,3);
            temp(i-nr,4) = Bdry(i+1,4);
            temp(i-nr,5) = 1;
            temp(i-nr,6) = 1;
            temp(i-nr+1,:) = []; % remove double entries
            nr = nr+1;
        end
    end
    
    Bdry = temp;
    
    N = length(Bdry(:,1));
    count = 0;
    
    % assembly of the tracematrix
    for i = 1:N_Omega
        for j = 1:N_Gamma
            
            val = 0;
            
            if ismember(i,position_array) % B_ij is only nonzero if we have a boundary node
                
                % get position of the i-th Omega node
                k1 = 1; % position of the node
                while ~(i == Bdry(k1,3))
                    k1 = k1+1; % search for the Node in bdry;
                end
                
                % get position of the j-th Gamma node
                k2 = 1; % position of the node
                while ~(j == Bdry(k2,4))
                    k2 = k2+1; % search for the Node in bdry;
                end
                
                if k1==k2 % we have the same nodes
                    p = 0;
                    n = 0;
                    
                    if k1-1 == 0
                        p = N;
                    else
                        p = k1-1;
                    end
                    
                    if k1+1 == N+1
                        n = 1;
                    else
                        n = k1+1;
                    end
                    
                    val = (norm(Bdry(p,1:2)-Bdry(k1,1:2))+norm(Bdry(n,1:2)-Bdry(k1,1:2)))/2; % the value of the next nodes is zero
                else % nodes are different; we need to get special intervalls
                
                    % get intervall of Omega-Node
                    p1 = k1-1; % position of prevous node
                
                    if p1 == 0
                        p1 = N; % don't get out of array
                    end
                
                    while Bdry(p1,5) == 0 % Gamma origin
                        p1 = p1-1;
                        if p1 == 0
                            p1 = N; % don't get out of array
                        end
                    end
                
                    n1 = k1+1; % position of next node
                    if n1 == N+1
                        n1 = 1; % don't get out of array
                    end
                
                    while Bdry(n1,5) == 0 % Gamma origin
                        n1 = n1+1;
                        if n1 == N+1
                            n1 = 1; % don't get out of array
                        end
                    end
                    
                    % get intervall of Gamma-Node
                    p2 = k2-1; % position of prevous node
                
                    if p2 == 0
                        p2 = N; % don't get out of array
                    end
                
                    while Bdry(p2,6) == 0 % Omega origin
                        p2 = p2-1;
                        if p2 == 0
                            p2 = N; % don't get out of array
                        end
                    end
                
                    n2 = k2+1; % position of next node
                    if n2 == N+1
                        n2 = 1; % don't get out of array
                    end
                
                    while Bdry(n2,6) == 0 % Omega origin
                        n2 = n2+1;
                        if n2 == N+1
                            n2 = 1; % don't get out of array
                        end
                    end
                    
                    % We now have intervalls [p1,n1] and [p2,n2] which are
                    % the supports of B\phi and \psi. We now need to
                    % intersect these and check if it is empty or not.
                    % If not then the supports overlap and we can calculate
                    % an integral over the intersection.
                    
                    % produce intervall arrays
                    
                    v1 = p1;
                    temp = p1;
                    while ~(temp == n1)
                        temp = temp+1;
                        if temp == N+1
                            temp = 1;
                        end
                        v1 = [v1,temp];
                    end
                    
                    v2 = p2;
                    temp = p2;
                    while ~(temp == n2)
                        temp = temp+1;
                        if temp == N+1
                            temp = 1;
                        end
                        v2 = [v2,temp];
                    end
                    
                    v = intersect(v1,v2);
                    
                    if length(v) == 2 % two points in the intersection
                        if ismember(p1,v) % then n2 is also in v
                            x1 = Bdry(p1,1:2);
                            x2 = Bdry(n2,1:2);
                            x = (x1+x2)/2; % choose center to evaluate
                                
                            val = norm(x1-x2) * norm(x-x1)/norm(Bdry(k1,1:2)-x1) * norm(x-x2)/norm(Bdry(k2,1:2)-x2); % approximate integral in the center because on the edges it is 0
                        else % the other way round (p2,n1)
                            x1 = Bdry(n1,1:2);
                            x2 = Bdry(p2,1:2);
                            x = (x1+x2)/2; % choose center to evaluate
                                
                            val = norm(x1-x2) * norm(x-x1)/norm(Bdry(k1,1:2)-x1) * norm(x-x2)/norm(Bdry(k2,1:2)-x2); % approximate integral in the center because on the edges it is 0
                        end
                    elseif length(v) == 3 % three points in the intersection
                        if ismember(k1,v) % then k2 is not in the intersection
                            x = Bdry(k2,1:2);
                            x1 = Bdry(p2,1:2);
                            x2 = Bdry(n2,1:2);
                            if norm(Bdry(k1,1:2)-x)/norm(x1-x) < norm(Bdry(k1,1:2)-x)/norm(x2-x)
                                temp = norm(Bdry(k1,1:2)-x)/norm(x1-x);
                            else
                                temp = norm(Bdry(k1,1:2)-x)/norm(x2-x);
                            end
                            val = temp * (norm(Bdry(v(1),1:2)-Bdry(k1,1:2))+norm(Bdry(v(3),1:2)-Bdry(k1,1:2)))/2;
                        elseif ismember(k2,v) % then k1 is not in the intersection
                            x = Bdry(k1,1:2);
                            x1 = Bdry(p1,1:2);
                            x2 = Bdry(n1,1:2);
                            if norm(Bdry(k2,1:2)-x)/norm(x1-x) < norm(Bdry(k2,1:2)-x)/norm(x2-x)
                                temp = norm(Bdry(k2,1:2)-x)/norm(x2-x);
                            else
                                temp = norm(Bdry(k2,1:2)-x)/norm(x2-x);
                            end
                            val = temp * (norm(Bdry(v(1),1:2)-Bdry(k2,1:2))+norm(Bdry(v(3),1:2)-Bdry(k2,1:2)))/2;
                        end
                    elseif length(v) == 4 % four points in the intersection + k1 and k2
                        x = Bdry(k2,1:2);
                        x1 = Bdry(p2,1:2);
                        x2 = Bdry(n2,1:2);
                        temp1 = 0;
                        if norm(Bdry(k1,1:2)-x)/norm(x1-x) < norm(Bdry(k1,1:2)-x)/norm(x2-x)
                            temp1 = norm(Bdry(k1,1:2)-x)/norm(x1-x);
                        else
                            temp1 = norm(Bdry(k1,1:2)-x)/norm(x2-x);
                        end
                        x = Bdry(k1,1:2);
                        x1 = Bdry(p1,1:2);
                        x2 = Bdry(n1,1:2);
                        temp2 = 0;
                        if norm(Bdry(k2,1:2)-x)/norm(x1-x) < norm(Bdry(k2,1:2)-x)/norm(x2-x)
                            temp2 = norm(Bdry(k2,1:2)-x)/norm(x1-x);
                        else
                            temp2 = norm(Bdry(k2,1:2)-x)/norm(x2-x);
                        end
                        
                        val = (temp1+temp2)/3 * (norm(Bdry(v(1),1:2)-Bdry(v(2),1:2)) + norm(Bdry(v(3),1:2)-Bdry(v(2),1:2)) + norm(Bdry(v(3),1:2)-Bdry(v(4),1:2)));
                    end
                end
            end
            B(i,j) = val;
        end
    end
    
end