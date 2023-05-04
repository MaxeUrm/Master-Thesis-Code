% Bulk potential derivative

% Concav part
function sol = W_minus_bulk(x)
    % 0-function
    %
    n = length(x);
    sol = zeros(n,1);
    %}
    
    % double-well
    %
    sol = - x;
    %}
end