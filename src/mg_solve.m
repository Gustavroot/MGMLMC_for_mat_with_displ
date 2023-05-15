function [x] = mg_solve(b,mgh,level_nr,tol,from_pgmres)

  % some MG params
  maxiter = 1;

  x = zeros(size(mgh.D{level_nr},1),1);
  b_norm = norm(b);

  for i=1:maxiter
    % compute residual
    r = b - mgh.D{level_nr}*x;
    rel_res = norm(r)/b_norm;
    if rel_res<tol
        break;
    end
    % get an approximation to the error
    mgh.bs{level_nr} = r(:);
    [mgh] = two_level_mg(mgh,level_nr,from_pgmres);
    % correct the solution
    x = x + mgh.xs{level_nr};
  end

  %iters = i;
end