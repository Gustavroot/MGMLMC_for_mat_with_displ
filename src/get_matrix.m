function [D] = get_matrix(filename)

  if filename~="matrices/D16/Dreord"
    error("Restricting for now to the example matrix in matrices/D16/Dreord");
  end

  fprintf("Getting matrix ...\n");
  t_start = tic;
  load(strcat(filename,"_top.mat"));
  load(strcat(filename,"_mid.mat"));
  load(strcat(filename,"_low.mat"));

  n = size(Dreord_top,2);
  D = sparse(n,n);
  D(1:n/3,1:n) = Dreord_top;
  D(n/3+1:2*n/3,1:n) = Dreord_mid;
  D(2*n/3+1:n,1:n) = Dreord_low;

  D = D - 0.01*speye(n);

  %tic;
  %d = eigs(D,1,'smallestreal','SubspaceDimension',50,'MaxIterations',50,'Tolerance',0.1);
  %toc
  %d

  fprintf("...done\n");
  t_end = toc(t_start);
  fprintf("Elapsed time : %f\n\n",t_end);

end