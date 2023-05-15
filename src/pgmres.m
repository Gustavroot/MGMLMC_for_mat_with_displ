function [x,iters] = pgmres(b,mgh,level_nr,tol)

  maxiter = 1000;

  % create anonymous function handle for preconditioner
  from_pgmres = 1;
  prec = @(bx) mg_solve(bx,mgh,level_nr,tol,from_pgmres);

  restart = 20;
  [x,~,~,iter] = gmres(mgh.D{level_nr},b,restart,tol,maxiter/restart,prec);

  iters = iter(1)*restart+iter(2);
end