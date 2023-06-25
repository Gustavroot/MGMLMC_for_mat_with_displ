% mgh : multigrid hierarchy
% alg_type : use either Hutchinson or MGMLMC
% k : number of deflation vectors
function [tracex,variance,iters] = compute_trace(k,mgh,alg_type,tol,maxiter,level_nr)

  % this function is used in case MGMLMC is chosen
  function [ex] = B(cx)
    ex = cx;
    for ix=(level_nr-1):-1:1
        ex = mgh.P{ix}*ex;
    end
    ex = mgh.GPM{1}'*(mgh.Ptilde{1}*(mgh.GPM{1}*ex));
    for ix=1:(level_nr-1)
        ex = mgh.R{ix}*ex;
    end
  end

  fprintf("Computing trace ...\n");
  tstart = tic;

  solver_tol = 1.0e-8;

  if alg_type=="Hutch"
    % handle for the operator to pass to Hutchinson
    A = @(bx) pgmres((mgh.GPM{1}'*(mgh.Ptilde{1}*(mgh.GPM{1}*bx))),mgh,1,solver_tol);
    % compute the variance
    if k>0
      [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
    else 
      [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{1},1));
    end
  elseif alg_type=="mgmlmc"
    % handle for the operator to pass to Hutchinson
    %B = @(cx) (mgh.P{1}'*(mgh.GPM{1}'*(mgh.Ptilde{1}*(mgh.GPM{1}*(mgh.P{1}*cx)))));
    % FIXME : double work being done here when computing E(bx)
    if level_nr==length(mgh.D)-1
      A = @(bx) ( pgmres(B(bx),mgh,level_nr,solver_tol) - ...
                  mgh.P{level_nr}*mgh.invD{level_nr+1}*(mgh.R{level_nr}*B(bx)));
    else
      A = @(bx) ( pgmres(B(bx),mgh,level_nr,solver_tol) - ...
                  mgh.P{level_nr}*pgmres(mgh.R{level_nr}*B(bx),mgh,level_nr+1,solver_tol));
    end
    % compute the variance and trace
    if k>0
      [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
    else 
      [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{2},1));
    end
  end
  tend = toc(tstart);
  fprintf("... done (elapsed time = %f)\n",tend);
end
