%At the first level we color the lattice matrix, so we have to extend the
%probing vectors to the Dirac matrix
function P = extend_probing_vectors(V, dof)

P = sparse(dof*size(V,1), size(V,2));

for i = 1:size(V,2)
    indices = find(V(:,i) == 1);
    for j = 1:length(indices)
        idx_start = (indices(j)-1)*dof + 1;
        idx_end = idx_start + dof - 1;
        P(idx_start:idx_end,i) = 1;
    end
end

end