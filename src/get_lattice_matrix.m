%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%EVALUATE THE ORDERING OF THE LATTICE POINTS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A is a QCD matrix, a 12x12 block (3 colors x 4 spins) in A is associated
%to each lattice site -> we have to compute the Frobenius norm of every
%12x12 block of A so that we know, by looking at the nonzero entries, the
%interacting sites of the lattice

%We can color this matrix instead of the original Dirac matrix and then
%extend every color by a factor of 12
%Example: substitute every red with 12 reds and every black with 12 blacks

function F = get_lattice_matrix(A, dof)
    t_start = tic;
    m = size(A, 1);

    n = floor(m / dof);
    
    % Initialize F as a sparse matrix
    F = sparse(n, n);

    for i = 1:n
        for j = 1:n
            % Select the current block
            current_block = A((i - 1) * dof + 1 : i * dof, (j - 1) * dof + 1 : j * dof);
            
            % Compute the Frobenius norm of the current block
            fro_norm = norm(current_block, 'fro');
            
            % Store the norm only if it is not zero
            if fro_norm ~= 0
                F(i, j) = fro_norm;
            end
        end
    end

    % spy(F); % optional: visualize the sparse matrix
    t_end = toc(t_start);
    fprintf("Elapsed time : %f\n",t_end);
end