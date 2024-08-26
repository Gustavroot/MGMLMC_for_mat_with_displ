function y = A_for_svd(x,tflag,P3D,mgh,tol)

  %Ax  = @(bx) P3D*( mgh.GPM{1}'*( pgmres(mgh.GPM{1}*(P3D'*(mgh.W{1}*bx)),mgh,1,tol) ) );
  %AxH = @(bx) mgh.W{1}'*( P3D*( mgh.GPM{1}'*( mgh.g5{1}*( pgmres(mgh.g5{1}*(mgh.GPM{1}*(P3D'*bx)),mgh,1,tol) ) ) ) );

  [~,G5_3D,~] = get_3D_params(mgh);

  Ax  = @(bx) A_hutch(bx,P3D,mgh.W{1},mgh,tol,"3D");
  AxH = @(bx) mgh.W{1}'*( G5_3D*( A_hutch( bx,P3D,G5_3D,mgh,tol,"3D" ) ) );

  if strcmp(tflag,'notransp')
    y = Ax(x);
  else
    y = AxH(x);
  end

end