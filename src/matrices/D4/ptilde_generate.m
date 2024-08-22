clc;clear;

load 4x4x4x4b6.0000id3n1.mat;

% invert D
Dinv = inv(D);

% build Ptilde, the matrix that displaces
nr_sites = 4*4*4*4;
dof = 12;
nr_displ_sites = 2*4;
Ptilde = sparse(nr_sites,nr_sites);
for j=1:nr_sites
  Ptilde(j,mod(j+nr_displ_sites-1,nr_sites)+1) = 1;
end
Ptilde = kron(Ptilde,speye(dof));

% extra term in variance
YY = ( Dinv' )*( Dinv.' );
XX = ( Dinv' * Ptilde )*( Dinv.' * (Ptilde').' );

%colormap('hot')
%imagesc(abs(XX))
%colorbar