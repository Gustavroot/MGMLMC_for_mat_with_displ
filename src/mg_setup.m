function [mgh] = mg_setup(D,nr_levels,nr_displ_sites)

  %global do_3D_traces;

  t_start = tic;

  mgh = struct("P",{{}},"R",{{}},"D",{{}},"g5",{{}}, ...
               "rs",{{}},"bs",{{}},"xs",{{}},"es",{{}}, ...
               "omega",{{}},"invD",{{}},"V",{{}},"Ptilde",{{}}, ...
               "GPM",{{}});

  fprintf("Constructing MG hierarchy ...\n")

  % dofs : [2,8,8,...]
  dofs = zeros(nr_levels);
  dofs(1) = 12;
  dofs(2) = 56;
  for i=3:nr_levels
    dofs(i) = 56;
  end
  mgh.dof = dofs;

  % aggrs : [4*4,2*2,2*2,...]
  aggrs = zeros(nr_levels-1);
  aggrs(1) = 4*4*4*4;
  for i=2:nr_levels-1
    aggrs(i) = 2*2*2*2;
  end
  mgh.aggr = aggrs;

  Df = D;
  mgh.D{1} = Df;
  for i=1:nr_levels-1

    nr_sites = size(mgh.D{i},1)/mgh.dof(i);
    mgh.g5{i} = kron(speye(nr_sites),blkdiag(speye(mgh.dof(i)/2),-speye(mgh.dof(i)/2)));

    [P,R,Dc] = build_2lvl_method_g5_comp(Df,dofs(i),dofs(i+1),aggrs(i),mgh.g5{i});

    mgh.P{i} = P;
    mgh.R{i} = R;
    mgh.D{i+1} = Dc;

    fprintf("\tDim of matrix at level %d : %d\n",i,size(mgh.D{i},1));
    fprintf("\tNumber of nonzero of matrix at level %d : %d\n",i,nnz(mgh.D{i}));

    % construction of Ptilde, which is the matrix giving the displacement on the lattice
    % IMPORTANT : this matrix permutes the chunks/sites on a vector upwards by an
    %             amount nr_displ_sites
    mgh.Ptilde{i} = sparse(nr_sites,nr_sites);
    for j=1:nr_sites
      mgh.Ptilde{i}(j,mod(j+nr_displ_sites-1,nr_sites)+1) = 1;
    end
    mgh.Ptilde{i} = kron(mgh.Ptilde{i},speye(mgh.dof(i)));

    %spy(mgh.Ptilde{i});
    %error("stop");

    if i==1
      load matrices/D16/GPM.mat;
      mgh.GPM{i} = GPM;
    end

    Df = Dc;
  end

  fprintf("\tDim of matrix at level %d : %d\n",nr_levels,size(mgh.D{nr_levels},1));
  fprintf("\tNumber of nonzero of matrix at level %d : %d\n",nr_levels,nnz(mgh.D{nr_levels}));
  nr_sites = size(mgh.D{nr_levels},1)/mgh.dof(nr_levels);
  mgh.g5{nr_levels} = kron(speye(nr_sites),blkdiag(speye(mgh.dof(nr_levels)/2),-speye(mgh.dof(nr_levels)/2)));

  % construction of Ptilde, which is the matrix giving the displacement on the lattice
  % IMPORTANT : this matrix permutes the chucks/sites on a vector upwards by an
  %             amount nr_displ_sites
  mgh.Ptilde{nr_levels} = sparse(nr_sites,nr_sites);
  for j=1:nr_sites
    mgh.Ptilde{nr_levels}(j,mod(j+nr_displ_sites-1,nr_sites)+1) = 1;
  end
  mgh.Ptilde{nr_levels} = kron(mgh.Ptilde{nr_levels},speye(mgh.dof(nr_levels)));

  % these buffers are needed in mg_solve.m
  for i=1:length(mgh.D)
    mgh.bs{i} = zeros(size(mgh.D{i},1),1);
    mgh.rs{i} = zeros(size(mgh.D{i},1),1);
    mgh.xs{i} = zeros(size(mgh.D{i},1),1);
    mgh.es{i} = zeros(size(mgh.D{i},1),1);
  end

  % compute parameters needed in Richardson, used as a smoother
  for i=1:length(mgh.D)-1
    if i==1
      lmax = 6.8436;
    else
      lmax = eigs(mgh.D{i},1,"largestabs",'Tolerance',1.0e-3);
      lmax = abs(lmax);
    end
    mgh.omega{i} = 1.0/(2.0*lmax/3.0);
  end

  % precompute the inverse at the coarsest level
  for i=1:length(mgh.D)
    if i<length(mgh.D)
      mgh.invD{i} = 0;
    else
      mgh.invD{i} = inv(mgh.D{i});
    end
  end

  fprintf("...done\n");

  t_end = toc(t_start);
  fprintf("Elapsed time : %f\n",t_end);

  fprintf("\n");

  % test correctness of MG hierarchy built
  mg_tests(mgh);

  fprintf("\n");
end
