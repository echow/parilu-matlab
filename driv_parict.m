clear all
% load ani3.mat % already diag scaled and rcm ordered
load /home/edmond/matrices.ufl/apache2.mat
a = Problem.A;
p = symrcm(a);
a = a(p,p);
a = diagscale(a);

n = length(a);
rng(0);
b = rand(n,1)-.5; 

numthreads = 4;
for numsweeps = 0:5

  u = triu(a);
  u = parict(a, u, 0, numsweeps, numthreads);

  [x flag relres iter resvec] = pcg(a, b, 1e-6, 2000, u', u);
  if flag, error('did not converge');
  else fprintf('%2d  %4d\n', numsweeps, iter); end

end


