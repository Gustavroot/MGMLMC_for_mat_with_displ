% trace_type : type of matrix for which the inverse is computed : inverse of D
%          or inverse of permuted D
% mgh : multigrid hierarchy
% alg_type : use either Hutchinson or MGMLMC
% k : number of deflation vectors
function [tracex,variance,iters] = compute_trace(trace_type,k,mgh,alg_type,tol,maxiter)

  fprintf("Computing trace ...\n");
  tstart = tic;

  % construct operator to pass to Block Power Iteration
  solver_tol = 1.0e-8;

  if alg_type=="Hutch"
    % with displaced lattice or not
    if trace_type=="invD*perm"
      fprintf("Computation of trace under construction.\n");
    end
    % handle for the operator to pass to Hutchinson
    A = @(bx) (mgh.GPM{1}'*(mgh.Ptilde{1}*(mgh.GPM{1}*pgmres(bx,mgh,1,solver_tol))));
    % compute the variance
    if k>0
      [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
    else 
      [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{1},1));
    end
  elseif alg_type=="mgmlmc"
    % with displaced lattice or not
    if trace_type=="invD*perm"
      fprintf("Computation of trace under construction.\n");
    end
    for i=1:length(mgh.D)-1
      % handle for the operator to pass to Hutchinson
      A = @(bx) ( pgmres(bx,mgh,1,solver_tol) - mgh.P{1}*pgmres(mgh.R{1}*bx,mgh,2,solver_tol) );
      % compute the variance
      if k>0
        [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
      else 
        [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{1},1));
      end
      % doing only the first level-difference for now
      break;
    end
    fprintf("Computation of trace under construction.\n");
  end
  tend = toc(tstart);
  fprintf("... done (elapsed time = %f)\n",tend);
end
