function y = A_mgmlmc_for_svd(x,tflag,P3D,mgh,solver_tol)

  % this function is only called for level_nr=1
  level_nr = 1;

  [~,G5_3D,~] = get_3D_params(mgh);

  Ax  = @(bx) A_mgmlmc(bx,P3D,mgh.W{1},mgh,solver_tol,level_nr);
  AxH = @(bx) mgh.W{1}'*( G5_3D*( A_mgmlmc(G5_3D*bx,P3D,mgh.W{1},mgh,solver_tol,level_nr) ) );

  if strcmp(tflag,'notransp')
    y = Ax(x);
  else

    y = AxH(x);
  end

end