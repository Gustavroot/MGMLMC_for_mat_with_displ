% cut matrices in half to be able to store on GitHub

% clc;clear;
% 
% load D.mat;
% n = size(D,1);
% D_top = D(1:n/3,1:n);
% D_mid = D(n/3+1:2*n/3,1:n);
% D_low = D(2*n/3+1:n,1:n);
% save('D_top.mat','D_top');
% save('D_mid.mat','D_mid');
% save('D_low.mat','D_low');

clc;clear;

load V.mat;
ni = size(V,1);
nj = size(V,2);
V_top = V(1:ni/4,1:nj);
V_mid1 = V(ni/4+1:2*ni/4,1:nj);
V_mid2 = V(2*ni/4+1:3*ni/4,1:nj);
V_low = V(3*ni/4+1:ni,1:nj);
save('V_top.mat','V_top');
save('V_mid1.mat','V_mid1');
save('V_mid2.mat','V_mid2');
save('V_low.mat','V_low');

% clc;clear;
% 
% load Dreord.mat;
% n = size(Dreord,1);
% Dreord_top = Dreord(1:n/3,1:n);
% Dreord_mid = Dreord(n/3+1:2*n/3,1:n);
% Dreord_low = Dreord(2*n/3+1:n,1:n);
% save('Dreord_top.mat','Dreord_top');
% save('Dreord_mid.mat','Dreord_mid');
% save('Dreord_low.mat','Dreord_low');