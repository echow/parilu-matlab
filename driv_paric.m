clear all
% load uni7
% p = sum(spones(a));
% j = (p ~= 1);
% a = a(j,j);
% p = symrcm(a);
% a = a(p,p);
% a = diagscale(a);

load /home/edmond/matrices.uniscaled/uni7scaled.mat

n = length(a);
%b = rand(n,1)-.5; 
b = load('b.txt');

%u = ichol(a)';
%u = triu(a);
%u = speye(n);

numthreads = 4;
for trial = 1:5
for numsweeps = 0:6

% pat = spones(iluk(a,2)');
% u = a.*pat + eps*pat;
  u = triu(a);

  u = paric(a, u, 1, numsweeps, numthreads);
  
  [x flag relres iter resvec] = pcg(a, b, 1e-6, 2000, u', u);
  if flag, error('did not converge'); 
  else fprintf('%2d  %4d\n', numsweeps, iter); end

end
end
