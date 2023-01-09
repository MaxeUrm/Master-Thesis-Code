% Approaching integrals for mass preservation and energy dissipation in 1d

function sol = integral_approx1d (val, Nodes, Elements)
    N = length(Elements);
    sol = 0;
    for i = 1:N
        
        % get baselength
        x1 = Nodes(Elements(i,1),:);
        x2 = Nodes(Elements(i,2),:);
        
        A_base = norm(x1-x2);
        
        % get mean heigth
        h = (val(Elements(i,1),1) + val(Elements(i,2),1))/2;
        
        % Integral over the element
        sol = sol + A_base * h;
    end
end