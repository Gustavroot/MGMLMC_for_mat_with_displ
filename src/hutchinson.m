function [tracex,variance,iters] = hutchinson(A,V,k,tol,maxiter,n,colors,level_nr,d,alg_type,mgh)

  ests = zeros(maxiter,1);
  variances = zeros(maxiter,1);
  fprintf("(stochs solves) ");
  
  if d > 0
    P = get_probing_vectors(colors);
    num_colors = size(P,2);
  
    if level_nr == 1
        dof = 12;
    else
        dof = 56;
    end
  
    P = extend_probing_vectors(P, dof);
  else
    num_colors = 1;
  end
  
  tracex = 0.0;
  variance = 0.0;

  for j = 1:num_colors
    for i=1:maxiter
        fprintf(".");
        
        if d > 0
            % Stochastic probing vectors
            z = build_spv(P(:,j)); 
        else
            % Rademacher vectors
            z = 2*randi([0 1],n,1)-1;
        end
   

        %MG-Def does deflation from the right
        if alg_type=="MG-Def"
            if k>0
                ests(i) = z'*(A(z)-(mgh.P{1}*(mgh.g5{2}*(mgh.U_k*(mgh.Lambda_hat_inv*(mgh.U_k'*(mgh.R{1}*z)))))));
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

        %fprintf("latest variance = %f\n",variances(i));
    end
    tracex = tracex + est_trace;
    variance = variance + variances(i);
  end
  fprintf("\n");

  iters = i;

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

  tracex = tracex + small_contr;
end

