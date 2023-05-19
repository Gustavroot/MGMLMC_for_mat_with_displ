clc; clear;

LASTN = maxNumCompThreads(1);

% some assumptions all over the code:
% -- tol for solves is always 1.0e-8
% -- the W-cycle and smoothing iters are set in two_level_mg.m
% -- when singular vectors are mentioned, it's of D and not D^{-1}

addpath("src/");

filename = "matrices/D16/Dreord";
D = get_matrix(filename);

% the shift on the lattice .. as we're doing here a lattice with a
% displacement. The displacement is on the second direction, that's why
% the 16 extra factor
nr_displ_sites = 5*(16^3);
%nr_displ_sites = 0;

% multigrid hierarchy
nr_levels = 3;
mgh = mg_setup(D,nr_levels,nr_displ_sites);

% compute and save singular values

Bxx = full(mgh.P{1}'*(mgh.GPM{1}'*mgh.Ptilde{1}*mgh.GPM{1})*mgh.P{1});
invD2 = inv(full(mgh.D{2}));
invD3 = inv(full(mgh.D{3}));
Cxx = Bxx*invD2;
Dxx = invD2*Bxx;
Exx = Bxx*( invD2 - full(mgh.P{2}*invD3*mgh.P{2}') );
Fxx = ( invD2 - full(mgh.P{2}*invD3*mgh.P{2}') )*Bxx;

fprintf("computing eig #1\n");
d1 = eig( Fxx'*Fxx );
fprintf("saving #1\n");
save('d1.mat','d1');

fprintf("computing eig #2\n");
d2 = eig( Exx'*Exx );
fprintf("saving #2\n");
save('d2.mat','d2');

fprintf("computing eig #3\n");
d3 = eig( Dxx'*Dxx );
fprintf("saving #3\n");
save('d3.mat','d3');

fprintf("computing eig #4\n");
d4 = eig( Cxx'*Cxx );
fprintf("saving #4\n");
save('d4.mat','d4');

fprintf("computing eig #5\n");
d5 = eig( invD2'*invD2 );
fprintf("saving #5\n");
save('d5.mat','d5');