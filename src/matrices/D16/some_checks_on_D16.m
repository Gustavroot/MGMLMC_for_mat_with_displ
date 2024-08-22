clc;clear;

LASTN = maxNumCompThreads(1);

load D16.mat;

nr_sites = size(D,1)/12;

for i=1:nr_sites
  jbeg = 1 + (i-1)*12;
  jend = jbeg + 11;
  norm_loc = norm(D(1:12,jbeg:jend),'fro');
  if norm_loc>0
    fprintf("id = %d\n",i);
  end
end
