% This function computes approximate deflation vectors via block power
% iteration
% param1 : type of deflation vectors : eigenvectors, right singular vectors
%          or left singular vectors
% param2 : number of deflation vectors
% param3 : multigrid hierarchy
% param4 : use either Hutchinson or MGMLMC
% param5 : iterations within Block Power Iteration
function [mgh] = compute_deflation_vectors(defl_type,k,mgh,alg_type,bpi_iters)

  fprintf("Constructing deflation vectors via BPI ...\n");
  tstart = tic;

  % construct operator to pass to Block Power Iteration
  tol = 1.0e-4;

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
    fprintf("Construction of deflation vectors under construction.\n");
  end
  tend = toc(tstart);
  fprintf("... done (elapsed time = %f)\n",tend);
end