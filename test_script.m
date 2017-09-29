numsweeps  = 50; % useful cases: 1 and small values and large values
numthreads = 4; % useful cases: 1 and 4 and larger values

a = delsq(numgrid('S',10));
a = diagscale(a);
u_ichol = ichol(a)';
[l_ilu u_ilu] = ilu(a);

% test paric_ref
u = triu(a);
u = paric_ref(a, u, numsweeps);
fprintf('paric_ref:     %e\n', norm(u-u_ichol,'fro')/norm(u_ichol,'fro'));

% test paric_ref with residual history computation
u = triu(a);
[u resid] = paric_ref(a, u, numsweeps);
fprintf('paric_ref:     %e\n', norm(u-u_ichol,'fro')/norm(u_ichol,'fro'));

% test paric/synchronous
u = triu(a);
u = paric(a, u, 1, numsweeps, numthreads);
fprintf('paric/sync:    %e\n', norm(u-u_ichol,'fro')/norm(u_ichol,'fro'));

% test paric/asynchronous
u = triu(a);
u = paric(a, u, 0, numsweeps, numthreads);
fprintf('paric/async:   %e\n', norm(u-u_ichol,'fro')/norm(u_ichol,'fro'));


% test parilu_ref
l = tril(a);
u = triu(a);
[l u] = parilu_ref(a, l, u, numsweeps);
fprintf('parilu_ref:L   %e\n', norm(l-l_ilu,'fro')/norm(l_ilu,'fro'));
fprintf('parilu_ref:U   %e\n', norm(u-u_ilu,'fro')/norm(u_ilu,'fro'));

% test parilu
l = tril(a);
u = triu(a);
[l u] = parilu(a, l, u, 0, numsweeps, numthreads);
fprintf('parilu/asy:L   %e\n', norm(l-l_ilu,'fro')/norm(l_ilu,'fro'));
fprintf('parilu/asy:U   %e\n', norm(u-u_ilu,'fro')/norm(u_ilu,'fro'));

% test pariutdt_ref (sync)
u = triu(a);
d = speye(length(a));
[u d] = pariutdu_ref(a, u, d, numsweeps);
fprintf('pariutdu_ref:  %e\n', norm(u-l_ilu','fro')/norm(l_ilu,'fro'));
