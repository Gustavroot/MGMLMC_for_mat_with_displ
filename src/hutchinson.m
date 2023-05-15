function [tracex,variance,iters] = hutchinson(A,V,k,tol,maxiter,n)

  ests = zeros(maxiter,1);
  variances = zeros(maxiter);
  fprintf("(stochs solves) ");
  for i=1:maxiter
    fprintf(".");
    z = rand(n,1);
    % deflation is from the left, as we're assuming using RSVs
    if k>0
      zdefl = z - V*(V'*z);
    else
      zdefl = z;
    end
    ests(i) = zdefl'*(A(z));

    % estimation of the trace
    est_trace = mean(ests(1:i));

    % estimation of the variance
    var_buff = ests(1:i)-est_trace;
    var_buff = abs(var_buff);
    var_buff = var_buff.^2;
    variances(i) = sum(var_buff)/i;

    fprintf("latest variance = %f\n",variances(i));

    % TODO : implement stop based on tolerance
  end
  fprintf("\n");

  % when using deflation, compute the trace of the small matrix
  if k>0
    W = zeros(n,k);
    fprintf("(direct solves) ");
    for i=1:k
      fprintf(".");
      W(:,i) = A(V(:,i));
    end
    fprintf("\n");
    small_contr = trace(V'*W);
  else
    small_contr = 0;
  end

  iters = i;
  tracex = est_trace + small_contr;
  variance = variances(i);
end
