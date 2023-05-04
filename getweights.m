% An Integration-function for Newton Cotes for second order konvex part

function weights = getweights(a,b,n)

    % we have to solve a linear equation
    
    h = (b-a)/n; % distance between equidistant nodes
    
    M = ones(n+1,n+1);
    
    for i = 1:(n+1)
        for j = 1:(n+1)
            M(i,j) = ((j-1)*h + a)^(i-1);
        end
    end
    
    v = zeros(n+1,1);
    
    for i = 1:(n+1)
        v(i) = 1/i * (b^i - a^i);
    end
    
    weights = M\v;
    
end