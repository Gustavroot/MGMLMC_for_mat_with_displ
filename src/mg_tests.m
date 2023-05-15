function [] = mg_tests(mgh)

  tests_ctr = 1;

  fprintf("Running tests on MG hierarchy and solver ...\n");

  % check if P is orthonormal at every level
  for i=1:length(mgh.P)
    rel_orth_norm = norm(mgh.P{i}'*mgh.P{i}-speye(size(mgh.P{i},2)),'fro')/norm(speye(size(mgh.P{i},2)),'fro');
    fprintf("Test #%d ... ",tests_ctr);
    try
      assert(rel_orth_norm<1.0e-15);
      fprintf("[PASSED]\n");
    catch
      fprintf("[FAILED] (rel error = %e)\n",rel_orth_norm);
    end
    tests_ctr = tests_ctr + 1;
  end

  % check if D is g5-Hermitian at every level
  for i=1:length(mgh.D)
    rel_g5Herm_norm = norm((mgh.g5{i}*mgh.D{i})'-(mgh.g5{i}*mgh.D{i}),'fro')/norm(mgh.g5{i}*mgh.D{i},'fro');
    fprintf("Test #%d ... ",tests_ctr);
    try
      assert(rel_g5Herm_norm<1.0e-15);
      fprintf("[PASSED]\n");
    catch
      fprintf("[FAILED] (rel error = %e)\n",rel_g5Herm_norm);
    end
    tests_ctr = tests_ctr + 1;
  end

  % testing lifted eigenvectors
  for i=1:length(mgh.D)-1
    Dc = mgh.D{i+1};
    [vc,~] = eigs(Dc,1,"smallestabs");
    vf = mgh.P{i}*vc;
    rq = (vf'*(mgh.D{i}*vf))/(vf'*vf);
    eig_res = norm(mgh.D{i}*vf-rq*vf)/norm(vf);
    fprintf("Test #%d ... ",tests_ctr);
    try
      assert(eig_res<1.0e-14);
      fprintf("[PASSED]\n");
    catch
      fprintf("[FAILED] (rel eig residual = %e)\n",eig_res);
    end
    tests_ctr = tests_ctr + 1;
  end

  % check if MG solver converges
  tol = 1.0e-8;
  for i=length(mgh.D)-1:-1:1
    b = rand(size(mgh.D{i},1),1);
    %[x,iters] = mg_solve(b,mgh,i,tol);
    tstart = tic;
    [x,iters] = pgmres(b,mgh,i,tol);
    tend = toc(tstart);
    rel_res_solve = norm(b-mgh.D{i}*x)/norm(b);
    fprintf("Test #%d ... ",tests_ctr);
    try
      assert(rel_res_solve<tol);
      fprintf("[PASSED], number of iterations = %d, time = %f\n",iters,tend);
    catch
      fprintf("[FAILED] (rel residual = %e, iters = %d, time = %f)\n",rel_res_solve,iters,tend);
    end
    tests_ctr = tests_ctr + 1;

  end

  fprintf("... done\n");

end