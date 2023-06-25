function [tracex,variance,iters] = hutchinson(A,V,k,tol,maxiter,n)

  ests = zeros(maxiter,1);
  variances = zeros(maxiter,1);
  fprintf("(stochs solves) ");
  for i=1:maxiter
    fprintf(".");

    % Rademacher vectors
    z = 2*randi([0 1],n,1)-1;

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

    %fprintf("latest variance = %f\n",variances(i));
  end
  fprintf("\n");

  iters = i;
  variance = variances(i);

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

  tracex = est_trace + small_contr;
end
