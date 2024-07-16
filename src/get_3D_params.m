function [P3D,G5_3D,rand_vec_size] = get_3D_params(mgh)

  % operator for projecting vectors from 4D to 3D
  dim4D = nthroot( size(mgh.D{1},1)/12,4 );
  size3D = size(mgh.D{1},1)/dim4D;
  size4D = size(mgh.D{1},1);
  P3D = sparse(size3D,size4D);
  % this variable, t, should be an input parameter
  t = 1;
  P3D( 1:size3D,1+(t-1)*size3D:size3D+(t-1)*size3D ) = speye(size3D);

  % size of random vectors in Hutchinson
  rand_vec_size = size(mgh.D{1},1)/dim4D;

  nr_sites_3D = rand_vec_size/12;
  G5_3D = kron(speye(nr_sites_3D),blkdiag(speye(12/2),-speye(12/2)));

end