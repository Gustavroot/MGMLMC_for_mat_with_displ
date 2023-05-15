function [P,R,Dc] = build_2lvl_method_g5_comp(Df,dof_f,dof_c,aggr_nr_sites)

  dof_aggr_f = dof_f*aggr_nr_sites;
  nr_aggrs = size(Df,1)/dof_aggr_f;
  nr_sites_per_aggr = dof_aggr_f/dof_f;

  ntv = dof_c/2;

  small_block = ones(dof_aggr_f,dof_c);

  Px = kron(speye(nr_aggrs),small_block);

  fprintf("\teigensolving ...\n");
  if size(Df,1)==(16*16*16*16*12)
    load matrices/D16/V_top.mat;
    load matrices/D16/V_mid.mat;
    load matrices/D16/V_low.mat;
    ni = size(V_top,1);
    V = ones(3*ni,ntv);
    V(1:ni,:) = V_top(:,1:ntv);
    V(ni+1:2*ni,:) = V_mid(:,1:ntv);
    V(2*ni+1:3*ni,:) = V_low(:,1:ntv);
  else
    [V,~,~] = eigs(Df,ntv,"smallestabs",'Tolerance',1.0e-4);
  end
  fprintf("\t... done\n");

  I_local = speye(dof_aggr_f);
  g5_local = ...
    kron(speye(nr_sites_per_aggr),blkdiag(speye(dof_f/2),-speye(dof_f/2)));

  % we could easily get rid of the following two for loops, via g5_local and an
  % extended g5

  for i=1:nr_aggrs

    % even

    ibeg = 1+(i-1)*dof_aggr_f;
    iend = i*dof_aggr_f;
    jbeg = 1 + (i-1)*dof_c;
    jend = jbeg+ntv-1;

    jbegx = 1;
    jendx = ntv;

    local_proj = (I_local-g5_local)/2;
    local_tv = V(ibeg:iend,jbegx:jendx);
    Px(ibeg:iend,jbeg:jend) = local_proj*local_tv;

    % odd

    ibeg = 1+(i-1)*dof_aggr_f;
    iend = i*dof_aggr_f;
    jbeg = 1 + (i-1)*dof_c+ntv;
    jend = jbeg+ntv-1;

    jbegx = 1;
    jendx = ntv;

    local_proj = (I_local+g5_local)/2;
    local_tv = V(ibeg:iend,jbegx:jendx);
    Px(ibeg:iend,jbeg:jend) = local_proj*local_tv;

  end

  P = kron(speye(nr_aggrs),small_block);

  for i=1:nr_aggrs
    ibeg = 1+(i-1)*dof_aggr_f;
    iend = i*dof_aggr_f;
    jbeg = 1 + (i-1)*dof_c;
    jend = i*dof_c;

    [Qx,~] = qr(full(Px(ibeg:iend,jbeg:jend)),0);
    P(ibeg:iend,jbeg:jend) = Qx;
  end

  R = P';
  Dc = R*(Df*P);

end
