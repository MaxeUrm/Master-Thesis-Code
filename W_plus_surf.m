% Surface potential derivative

% Convex part
function sol = W_plus_surf(x)
    % 0-function
    %{
    n = length(x);
    sol = zeros(n,1);
    %}
    
    % double-well
    %
    sol = x.^3;
    %}
end