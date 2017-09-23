function u = parict(a, u, numsweeps, numthreads)
% u = parict(a, u, numsweeps, numthreads)

% could redesign this function so that it only performs a single sweep

% notes:
% ordering for find is important
% 0-based indexing for mex file

n = length(a);
nz = nnz(triu(a));

for sweep = 1:numsweeps

  % find candidates (here, set c includes existing nonzeros)
  pat = spones(u);
  c = triu(spones(pat'*pat) + spones(a));

  % add all possible nonzeros to u
  u = u + eps*c;

  % run parilu step
  pat = spones(u); % updated pattern of u
  [Si Sj Sa] = find(a.*pat + eps*pat);
  u = paric_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, 1, numthreads);

  % delete small entries (keep only large entries)
  [ci cj cval] = find(u);
  [~, ind] = sort(abs(cval), 1, 'descend');
  u = sparse(ci(ind(1:nz)),cj(ind(1:nz)),cval(ind(1:nz)),n,n);

  % run parilu step
  pat = spones(u);
  [Si Sj Sa] = find(a.*pat + eps*pat);
  u = paric_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, 1, numthreads);

end

