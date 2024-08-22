function [tracex,variance,iters] = hutchinson(A,V,k,tol,maxiter,n, alg_type, mgh)

  ests = zeros(maxiter,1);
  variances = zeros(maxiter,1);
  fprintf("(stochs solves) ");
  for i=1:maxiter
    fprintf(".");

    % Rademacher vectors
    z = 2*randi([0 1],n,1)-1;

    %MG-Def does deflation from the right
    if alg_type=="MG-Def"
      if k>0
        ests(i) = z'*(A(z)-(mgh.P{1}*(mgh.g5{2}*...
             (mgh.U_k*(mgh.Lambda_hat_inv*(mgh.U_k'*(mgh.R{1}*z)))))));
      else
        ests(i) = z'*(A(z));
      end
      
    else
      % deflation is from the left, as we're assuming using RSVs
      if k>0
        zdefl = z - V*(V'*z);
      else
        zdefl = z;
      end
      ests(i) = zdefl'*(A(z));
    end
    % estimation of the trace
    est_trace = mean(ests(1:i));

    % estimation of the variance
    var_buff = ests(1:i)-est_trace;
    var_buff = abs(var_buff);
    var_buff = var_buff.^2;
    variances(i) = sum(var_buff)/i;

    fprintf("latest variance = %f\n",variances(i));
    fprintf("latest deflated trace = %f\n",real(est_trace));
  end
  fprintf("\n");

  iters = i;
  variance = variances(i);

  % when using deflation, compute the trace of the small matrix
  if k>0
    if alg_type=="MG-Def"
      small_contr = 0;
      for i=1:size(mgh.U_k,2)
        small_contr = small_contr + ...
                      mgh.U_k(:,i)'*mgh.g5{2}*(mgh.U_k(:,i)*mgh.Lambda_hat_inv(i,i));
      end
    else
      W = zeros(n,k);
      fprintf("(direct solves) ");
      for i=1:k
        fprintf(".");
        W(:,i) = A(V(:,i));
      end
      fprintf("\n");
      small_contr = trace(V'*W);
    end
  else
    small_contr = 0;
  end

  tracex = est_trace + small_contr;
end
