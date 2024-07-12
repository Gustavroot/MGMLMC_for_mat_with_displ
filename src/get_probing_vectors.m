% Build probing vectors for D given the coloring
%n is the dimension of D
function P = get_probing_vectors(colors)

num_colors = max(colors);
P = sparse(size(colors,1), num_colors);

for i = 1:num_colors
    indices = find(colors == i);
    P(indices, i) = 1;
end

end