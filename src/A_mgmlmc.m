% this function is used in case MGMLMC is chosen from within
% compute_trace.m, this is when doing 3D traces
function [ex] = A_mgmlmc(cx,P3D,mgh,solver_tol,level_nr)
  ex = cx;
  ex = P3D'*ex;
  ex = mgh.GPM{1}*ex;
  
  for ix=1:level_nr-1
    ex = mgh.R{ix}*ex;
  end
  
  if level_nr==length(mgh.D)-1
    ex = pgmres(ex,mgh,level_nr,solver_tol) - ...
         mgh.P{level_nr}*mgh.invD{level_nr+1}*(mgh.R{level_nr}*ex);
  else
    ex = pgmres(ex,mgh,level_nr,solver_tol) - ...
         mgh.P{level_nr}*pgmres(mgh.R{level_nr}*ex,mgh,level_nr+1,solver_tol);
  end
  
  for ix=level_nr-1:-1:1
    ex = mgh.P{ix}*ex;
  end
  
  ex = mgh.GPM{1}'*ex;
  ex = P3D*ex;
end