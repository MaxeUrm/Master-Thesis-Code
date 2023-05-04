% This function computes the maximal width of the Elements

function sol = max_width(Nodes_bulk, Elements_bulk, No_Elements)

    sol = 0;
    
    for i = 1:No_Elements
        x1 = Nodes_bulk(Elements_bulk(i,1),:);
        x2 = Nodes_bulk(Elements_bulk(i,2),:);
        x3 = Nodes_bulk(Elements_bulk(i,3),:);
        
        sol = max([sol; norm(x1-x2); norm(x1-x3); norm(x2-x3)]);
    end

end