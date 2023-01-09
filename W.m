% Potential function in the bulk

function sol = W(x)
    % 0-function
    %{
    n = length(x);
    sol = zeros(n,1);
    %}
    
    % double well potential
    %
    sol = 1/4 * (x.^4 - 2 * x.^2 + 1);
    %}

end