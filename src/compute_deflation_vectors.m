% This function computes approximate deflation vectors via block power
% iteration
% param1 : type of deflation vectors : eigenvectors, right singular vectors
%          or left singular vectors
% param2 : number of deflation vectors
% param3 : multigrid hierarchy
% param4 : use either Hutchinson or MGMLMC
% param5 : iterations within Block Power Iteration
function [mgh] = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters)

  global do_3D_traces;
  
  fprintf("Constructing deflation vectors via BPI ...\n");
  tstart = tic;

  tol = 1.0e-4;

  if do_3D_traces==1

    if defl_type~="RSVs"
      error("Only supporting deflation of singular vectors at the moment\n");
    end

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

    %solver_tol = 1.0e-8;

    nr_sites_3D = rand_vec_size/12;
    G5_3D = kron(speye(nr_sites_3D),blkdiag(speye(12/2),-speye(12/2)));

    % handle for the operator to pass to Hutchinson
    A = @(bx) G5_3D*( P3D*( mgh.GPM{1}'*( pgmres(mgh.GPM{1}*(P3D'*bx),mgh,1,tol) ) ) );

    mgh.V{1} = bpi(A,k,bpi_iters,rand_vec_size);

  else

    if alg_type=="Hutch"
      % the following matrix allows using the Hermitian operators
      if defl_type=="EVs"
        G = speye(size(D,1));
      else
        G = mgh.g5{1};
      end
      % handle for the operator to pass to BPI
      A = @(bx) pgmres(G*bx,mgh,1,tol);
      % call Block Power Iteration. In the case of not calling with "EVs",
      % this function returns approximate left singular vectors of D^{-1},
      % which are approximate right singular vectors of D
      mgh.V{1} = bpi(A,k,bpi_iters,size(mgh.D{1},1));
    elseif alg_type=="mgmlmc"
      for i=1:length(mgh.D)-1
        % the following matrix allows using the Hermitian operators
        if defl_type=="EVs"
          G = speye(size(mgh.D{1},1));
        else
          G = mgh.g5{1};
        end
        % handle for the operator to pass to BPI
        A = @(bx) ( pgmres(G*bx,mgh,1,tol) - mgh.P{1}*pgmres(mgh.R{1}*G*bx,mgh,2,tol) );
        % call Block Power Iteration. In the case of not calling with "EVs",
        % this function returns approximate left singular vectors of D^{-1},
        % which are approximate right singular vectors of D
        mgh.V{1} = bpi(A,k,bpi_iters,size(mgh.D{1},1));
  
        % only first level-difference for now
        break;
      end
      %fprintf("Construction of deflation vectors under construction.\n");
    end

  end

  tend = toc(tstart);
  fprintf("... done (elapsed time = %f)\n",tend);

end