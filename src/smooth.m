function [x] = smooth(D,b,iters,omega)

  % initial guess
  x = zeros(size(D,1),1);

  % Richardson
  for i=1:iters
    % compute the residual
    r = b - D*x;
    % update the solution
    x = x + omega*r;
  end

end