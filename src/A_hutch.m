function y = A_hutch(x,P3D,F,mgh,tol,dim)

  if strcmp(dim,'3D')
    y = P3D*( mgh.GPM{1}'*( pgmres(mgh.GPM{1}*(P3D'*(F*x)),mgh,1,tol) ) );
  else
    % this second case is when evaluating 4D traces
    error("Refactoring of 4D traces under construction\n");
  end

end