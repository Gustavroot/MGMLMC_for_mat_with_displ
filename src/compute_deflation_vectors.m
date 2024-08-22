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

  fprintf("Constructing deflation vectors ...\n");
  tstart = tic;

  % tolerance of the solves involved in all these calls to BPI
  tol = 1.0e-8;

  if do_3D_traces==1

    if defl_type~="RSVs"
      error("Only supporting deflation of singular vectors at the moment\n");
    end

    [P3D,G5_3D,~] = get_3D_params(mgh);

    % check if mgh.W{1} is unitary
    unit_norm = norm(mgh.W{1}'*mgh.W{1}-speye(size(mgh.W{1},1)),'fro')/...
                norm(speye(size(mgh.W{1},1)),'fro');

    % if mgh.W{1} is unitary, then we can use a G5_3D-Hermitian operator
    % for the SVD extraction
    if unit_norm<1.0e-15

      if alg_type=="Hutch"
        A_for_eig = @(bx) A_hutch(bx,P3D,G5_3D,mgh,tol,"3D");
        %A = @(bx) P3D*( mgh.GPM{1}'*( pgmres(mgh.GPM{1}*(P3D'*(G5_3D*bx)),mgh,1,tol) ) );
        %mgh.V{1} = bpi(A,k,bpi_iters,rand_vec_size);
      else
        A_for_eig  = @(bx) A_mgmlmc(bx,P3D,G5_3D,mgh,tol,1);
      end

      % in case mgh.W{1} is unitary, use eigs
      [X,~,flag] = eigs(A_for_eig,size(P3D,1),k,'largestabs', ...
                        'Tolerance',1.0e-5, ...
                        'MaxIterations',1000, ...
                        'SubspaceDimension',60 ...
                        );

    else

      if alg_type=="Hutch"
        A = @(x,tflag) A_for_svd(x,tflag,P3D,mgh,tol);
      else
        % this case is called only for the first difference level
        A = @(x,tflag) A_mgmlmc_for_svd(x,tflag,P3D,mgh,tol);
      end

      % in case mgh.W{1} is non-unitary, use svds
      [X,~,~,flag] = svds( @(x,tflag) A(x,tflag), ...
                           [size(P3D,1),size(P3D,1)], k, 'largest', ...
                           'Tolerance',1.0e-5, ...
                           'MaxIterations',1000, ...
                           'SubspaceDimension',60 ...
                           );
    end

    mgh.V{1} = X;
    if flag==0
      fprintf("all the requested singular vectors (for deflation) converged!\n");
    else
      error("not all the requested singular vectors (for deflation) converged!\n");
    end

  elseif do_3D_traces==0 && alg_type=="MG-Def"

    Q= full(mgh.g5{2}*mgh.D{2});
    [mgh.V{2}, mgh.Lambda{2}] = eigs(Q, k, 'smallestabs', 'Tolerance', 1e-2,'SubspaceDimension', 2*k, 'MaxIterations', 1000, 'Display', 1);
    mgh.U{2} = mgh.g5{2}*mgh.V{2}*sign(mgh.Lambda{2});  
    %Factorization to simplify calculations
    B = (mgh.U{2}'*mgh.D{2})*(mgh.g5{2}*mgh.U{2});
    [mgh.U_hat,mgh.Lambda_hat_inv] = eig(B);
    mgh.Lambda_hat_inv = inv(mgh.Lambda_hat_inv);
    mgh.U_k = mgh.U{2}* mgh.U_hat;
      
  else

    error("Refactoring of 4D deflation vectors under construction\n");

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
