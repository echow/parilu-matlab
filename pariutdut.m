function [u d] = pariutdut(a, u, d, sync_notsync, numsweeps, numthreads)
% [u d] = pariutdut(a, u, d, sync_notsync, numsweeps, numthreads)

% notes:
% ordering for find is important
% 0-based indexing for mex file
% could redesign this function so that it only performs a single sweep

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
  if (sync_notsync)
    [u d] = pariutdu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, 1, numthreads);
  else
    [u d] = pariutdu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, 1, numthreads);
  end

  % delete small entries (keep only large entries)
  [ci cj cval] = find(u);
  [~, ind] = sort(abs(cval), 1, 'descend');
  u = sparse(ci(ind(1:nz)),cj(ind(1:nz)),cval(ind(1:nz)),n,n);
  if nnz(diag(u) == 0)
    error('diagonal element dropped');
  end

  % run parilu step
  pat = spones(u);
  [Si Sj Sa] = find(a.*pat + eps*pat);
  if (sync_notsync)
    [u d] = pariutdu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, 1, numthreads);
  else
    [u d] = pariutdu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, 1, numthreads);
  end

end

