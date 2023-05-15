function [V] = bpi(A,k,s,n)

  % A : operator
  % k : number of vectors
  % s : number of Block Power Iteration steps
  % n : size of the system matrix

  V = rand(n,k);
  W = zeros(n,k);

  for i=1:s
    for j=1:k
      W(:,j) = A(V(:,j));
    end
    [V,~] = qr(W,0);
  end

end