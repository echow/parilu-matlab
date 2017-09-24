clear all
numthreads = 4;
first = 60;
for matnum=first:69
  filename = sprintf('/home/edmond/matrices.topopt/system_%05d.mat', matnum);
  prob = load(filename);
  a = prob.A;
  a = (a+a')/2;
  a = diagscale(a);
  afun = @(x) a'*x;
  rng(0);
  n = length(a);
  b = rand(length(a),1)-.5;
  fprintf('============ matrix %5d:  number of nonzeros %d\n', matnum, nnz(a));

  if matnum == first
    u = ichol(a)';
    continue;
  else
    for sweep = 0:10
      if (sweep > 0)
        u = parict(a, u, 1, numthreads);
      end

      [x flag relres iter resvec] = pcg(afun, b, 1e-6, 2000, u', u);
      if flag, error('did not converge');
      else fprintf('sweep %2d  %4d\n', sweep, iter); end
    end
  end
end

