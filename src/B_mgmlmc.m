% this function is used in case MGMLMC is chosen, for 4D traces
function [ex] = B_mgmlmc(cx,mgh,level_nr)
  ex = cx;
  for ix=(level_nr-1):-1:1
    ex = mgh.P{ix}*ex;
  end
  ex = mgh.GPM{1}*(mgh.Ptilde{1}*(mgh.GPM{1}'*ex));
  for ix=1:(level_nr-1)
    ex = mgh.R{ix}*ex;
  end
end