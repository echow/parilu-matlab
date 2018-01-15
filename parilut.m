function [l u] = parilut(a, l, u, sync_notsync, numsweeps, numthreads)
% [l u] = parilut(a, l, u, sync_notsync, numsweeps, numthreads)

% notes:
% ordering for find is important
% 0-based indexing for mex file
% could redesign this function so that it only performs a single sweep

n = length(a);
nzl = nnz(l); % bug fixed 11/23/2017;
nzu = nnz(u); % bug fixed 11/23/2017;

for sweep = 1:numsweeps

  % find candidates (here, set c includes existing nonzeros)
  c = spones( spones(l)*spones(u) + spones(a) );

  % add all possible nonzeros to l and u
  l = l + eps*tril(c);
  u = u + eps*triu(c);

  % run parilu step
  pat = spones(spones(l) + spones(u)); % updated pattern of l and u
  [Si Sj Sa] = find(a.*pat + eps*pat);
  l = l';
  if (sync_notsync)
    [l u] = parilu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, 1, numthreads);
  else
    [l u] = parilu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, 1, numthreads);
  end
  l = l';

  % delete small entries (keep only large entries)
  [ci cj cval] = find(l);
  [~, ind] = sort(abs(cval), 1, 'descend');
  l = sparse(ci(ind(1:nzl)),cj(ind(1:nzl)),cval(ind(1:nzl)),n,n);
  if nnz(diag(l) == 0)
    error('L diagonal element dropped');
  end

  [ci cj cval] = find(u);
  [~, ind] = sort(abs(cval), 1, 'descend');
  u = sparse(ci(ind(1:nzu)),cj(ind(1:nzu)),cval(ind(1:nzu)),n,n);
  if nnz(diag(u) == 0)
    error('U diagonal element dropped');
  end

  % run parilu step
  pat = spones(spones(l) + spones(u));
  [Si Sj Sa] = find(a.*pat + eps*pat);
  l = l';
  if (sync_notsync)
    [l u] = parilu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, 1, numthreads);
  else
    [l u] = parilu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, 1, numthreads);
  end
  l = l';

end

