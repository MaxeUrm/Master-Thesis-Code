% Bulk potential second derivative

% Convex part
function sol = W_plus_bulk_2(x)
    % 0-function
    %{
    n = length(x);
    sol = zeros(n,1);
    %}
    
    % double-well
    %
    sol = 3 * x.^2;
    %}
end