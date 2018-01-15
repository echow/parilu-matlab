function [u d] = pariutdu(a, u, d, sync_notsync, numsweeps, numthreads)
% [u d] = pariutdu(a, u, d, sync_nosync, numsweeps, numthreads)

% notes:
% ordering for find is column-major; important for async case
% 0-based indexing used for mex program

pat = spones(u);
[Si Sj Sa] = find(a.*pat + eps*pat);
if (sync_notsync)
  [u d] = pariutdu_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, numsweeps, numthreads);
else
  [u d] = pariutdu_async_mex(int32(Si)-1, int32(Sj)-1, Sa, u, d, numsweeps, numthreads);
end



