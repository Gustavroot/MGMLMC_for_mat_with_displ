function colors = graph_coloring(mgh, nr_levels, d)

assert(d > 0, "The coloring distance must be positive");

t_start = tic;
addpath("lattice_matrices/")
%create an array of array which at every level of the mg hiearchy contains
%the colors of the graph associated to the Dirac matrix at that level
colors = cell(nr_levels-1, 1);
fprintf("Graph coloring...\n");

%loop 1:nr_levels - 1 because the trace at the last level is computed
%exactly
for i = 1:nr_levels-1
    
    fprintf("Coloring the graph at the level number %d ...\n", i);
     
    if i <= 2 %because we already computed A = get_lattice_matrix(mgh.D{i},dof); for i = 1,2
        
        %At the first level we load the lattice matrix (which ignores the spin and
        %color dof) -> we color this matrix and the we extend the probing
        %vectors in Hutchinson.m to the Dirac matrix dimension 
        
        string = sprintf("level_%d_lattice_MGMLMC_16^4_Matrix.mat", i);
        load(string);
    else
        %for level_nr > 1 dof = 56
        A = get_lattice_matrix(mgh.D{i},56);
    end
    
    %We build the distance-d coloring of the graph of A by performing the
    %distance-1 coloring of A^d
    if d > 1
        A = sparseMatrixPower(A, d);
    end
    
    colors{i} = distance1Coloring(A);
        
    
    
    fprintf("Level number %d done...\n", i);
end

t_end = toc(t_start);
fprintf("... done (elapsed time = %f)\n",t_end);
end