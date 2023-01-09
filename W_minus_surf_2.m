% Surface potential second derivative

% Concav part
function sol = W_minus_surf_2(x)
    % 0-function
    %
    n = length(x);
    sol = zeros(n,1);
    %}
    
    % double-well
    %
    sol = - ones(n,1);
    %}
end