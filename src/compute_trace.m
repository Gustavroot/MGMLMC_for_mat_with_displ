% mgh : multigrid hierarchy
% alg_type : use either Hutchinson or MGMLMC
% k : number of deflation vectors
function [tracex,variance,iters] = compute_trace(k,mgh,alg_type,tol,maxiter,level_nr)

  global do_3D_traces;

  if do_3D_traces==1

    [P3D,~,rand_vec_size] = get_3D_params(mgh);

    fprintf("Computing trace ...\n");
    tstart = tic;

    solver_tol = 1.0e-8;

    if alg_type=="Hutch"
      % handle for the operator to pass to Hutchinson
      %A = @(bx) P3D*( mgh.GPM{1}'*( pgmres(mgh.GPM{1}*(P3D'*(mgh.W{1}*bx)),mgh,1,solver_tol) ) );
      A = @(bx) A_hutch(bx,P3D,mgh.W{1},mgh,tol,"3D");

      % compute the variance
      if k>0
        [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,rand_vec_size);
      else
        [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,rand_vec_size);
      end
    elseif alg_type=="mgmlmc"
      if level_nr==length(mgh.D)
        tracex = trace(A_mgmlmc_coarsest(0,P3D,mgh,0,level_nr));
        variance = 0.0;
        iters = 0;
      else
        A = @(bx) A_mgmlmc(bx,P3D,mgh.W{1},mgh,solver_tol,level_nr);

        % compute the variance and trace
        if k>0
          % for now, only the first difference level undergoes deflation
          if level_nr==1
            [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,rand_vec_size);
          else
            [tracex,variance,iters] = hutchinson(A,0,0,tol,maxiter,rand_vec_size);
          end
        else 
          [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,rand_vec_size);
        end
      end
    end
    tend = toc(tstart);
    fprintf("... done (elapsed time = %f)\n",tend);

  else

    error("Refactoring of 4D traces is under construction\n");
  
    B = @(cx) B_mgmlmc(cx,mgh,level_nr);

    fprintf("Computing trace ...\n");
    tstart = tic;

    solver_tol = 1.0e-8;

    if alg_type=="Hutch"
      % handle for the operator to pass to Hutchinson
      A = @(bx) pgmres((mgh.GPM{1}*(mgh.Ptilde{1}*(mgh.GPM{1}'*bx))),mgh,1,solver_tol);

      % compute the variance
      if k>0
        [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
      else
        [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{1},1));
      end
    elseif alg_type=="mgmlmc"
      if level_nr==length(mgh.D)
        tracex = trace(mgh.invD{level_nr});
        variance = 0.0;
        iters = 0;
      else
        % handle for the operator to pass to Hutchinson
        if level_nr==length(mgh.D)-1
          A = @(bx) ( pgmres(B(bx),mgh,level_nr,solver_tol) - ...
                      mgh.P{level_nr}*mgh.invD{level_nr+1}*(mgh.R{level_nr}*B(bx)));
        else
          A = @(bx) ( pgmres(B(bx),mgh,level_nr,solver_tol) - ...
                      mgh.P{level_nr}*pgmres(mgh.R{level_nr}*B(bx),mgh,level_nr+1,solver_tol));
        end
        % compute the variance and trace
        if k>0
          error("Deflation with MGMLMC is disabled at the moment\n")
          [tracex,variance,iters] = hutchinson(A,mgh.V{1},k,tol,maxiter,size(mgh.D{1},1));
        else 
          [tracex,variance,iters] = hutchinson(A,0,k,tol,maxiter,size(mgh.D{level_nr},1));
        end
      end
    end
    tend = toc(tstart);
    fprintf("... done (elapsed time = %f)\n",tend);
  end

end
