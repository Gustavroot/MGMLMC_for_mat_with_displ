% this function is used in case MGMLMC is chosen from within
% compute_trace.m, this is when doing 3D traces
function [M] = A_mgmlmc_coarsest(cx,P3D,mgh,solver_tol,level_nr)

  % important : this function is called for 3D traces only at the moment

  X1 = mgh.W{1};
  X1 = P3D'*X1;
  X1 = mgh.GPM{1}*X1;

  for ix=1:level_nr-1
    X1 = mgh.R{ix}*X1;
  end

  %X1 = mgh.invD{level_nr}*X1;

  X2 = mgh.P{level_nr-1};
  for ix=level_nr-2:-1:1
    X2 = mgh.P{ix}*X2;
  end

  X2 = mgh.GPM{1}'*X2;
  X2 = P3D*X2;

  M = mgh.invD{level_nr}*(X1*X2);
end