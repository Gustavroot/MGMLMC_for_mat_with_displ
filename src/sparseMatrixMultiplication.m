function C = sparseMatrixMultiplication(A, B)
    % Ensure A and B are sparse matrices
    if ~issparse(A)
        A = sparse(A);
    end
    if ~issparse(B)
        B = sparse(B);
    end
    % Perform sparse matrix multiplication
    C = A * B;
end
