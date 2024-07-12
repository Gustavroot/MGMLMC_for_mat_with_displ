function col = distance1Coloring(A)
    % Check if the input matrix A is square
    [n, m] = size(A);
    if n ~= m
        error('The adjacency matrix must be square.');
    end
    
    t_start = tic;
    
    fprintf("Coloring the graph...\n");
    %Initialize the nodes array
    K = (1:n);
    
    % Initialize the color array
    col = zeros(n, 1);
    
    % Iterate over the sequence of nodes K
    for i = 1:n
        %fprintf("Coloring node %d\n", i);
        node = K(i);
        
        % Find the colors of the adjacent nodes
        neighborColors = col(A(node, :) > 0);
        
        % Find the minimum color that is not used by the neighbors
        availableColor = 1;
        while any(neighborColors == availableColor)
            availableColor = availableColor + 1;
        end
        
        % Assign the color to the current node
        col(node) = availableColor;
    end
    
    fprintf("...done\n");
    t_end = toc(t_start);
    fprintf("Elapsed time : %f\n",t_end);
end