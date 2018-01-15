function [l u] = parilu(a, l, u, sync_notsync, numsweeps, numthreads)
% [l u] = parilu(a, l, u, sync_notsync, numsweeps, numthreads)

pat = spones(spones(l) + spones(u));
[Si Sj Sa] = find(a.*pat + eps*pat);

l = l'; % convert l to row-based
if (sync_notsync)
  [l u] = parilu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, numsweeps, numthreads);
else
  [l u] = parilu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, numsweeps, numthreads);
end
l = l'; % convert back


