clc;clear;
LASTN = maxNumCompThreads(1);

load D16/D16_reord.mat;

% let's build a block-diagonal preconditioner

tic;
fprintf("Constructing preconditioner ...\n");
chnk = 48;
B = kron(speye(size(Dreord,1)/chnk),ones(chnk));
for i=1:size(Dreord,1)/chnk
  if mod(i,1000)==0
    fprintf(" i = %d ... ",i);
  end
  idx = 1 + (i-1)*chnk;
  B(idx:idx+(chnk-1),idx:idx+(chnk-1)) = inv(Dreord(idx:idx+(chnk-1),idx:idx+(chnk-1)));
end
fprintf("\n");
Bfun = @(bx) B*bx;
fprintf("... done\n");
toc

Afun = @(bx) gmres(Dreord,bx,40,1.0e-12,100,Bfun);

ntv = 28;
tol = 1.0e-4;

tic;
[V,d] = eigs(Afun,size(Dreord,1),ntv,"smallestabs",'Tolerance',tol,'SubspaceDimension',60,...
               'MaxIterations',500,'Display',1);
toc

save('D16/Vy.mat','V');
