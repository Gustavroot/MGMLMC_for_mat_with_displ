function Apow = sparseMatrixPower(A, n)
    % Ensure A is a sparse matrix
    if ~issparse(A)
        A = sparse(A);
    end
    
    % Initialize Apow as the identity matrix of the same size as A
    Apow = speye(size(A));
    
    % Perform n-th power using sparse matrix multiplication
    for i = 1:n
        Apow = sparseMatrixMultiplication(Apow, A);
    end
end