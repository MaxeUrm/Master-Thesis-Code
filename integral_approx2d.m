% Approaching integrals for mass preservation and energy dissipation in 2d

function sol = integral_approx2d (val, Nodes, Elements)
    N = length(Elements);
    sol = 0;
    for i = 1:N
        
        % get basearea
        x1 = Nodes(Elements(i,1),:);
        x2 = Nodes(Elements(i,2),:);
        x3 = Nodes(Elements(i,3),:);
        
        M = [x1;x2;x3];
        
        A_base = polyarea(M(:,1),M(:,2));
        
        % get mean heigth
        h = (val(Elements(i,1),1) + val(Elements(i,2),1) + val(Elements(i,3),1))/3;
        
        % Integral over the element
        sol = sol + A_base * h;
    end
end