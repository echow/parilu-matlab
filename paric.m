function u = paric(a, u, sync_notsync, numsweeps, numthreads)
% u = paric(a, u, sync_nosync, numsweeps, numthreads)

% notes:
% ordering for find is column-major; important for async case
% 0-based indexing used for mex program

pat = spones(u);
[Si Sj Sa] = find(a.*pat + eps*pat);
if (sync_notsync)
  u = paric_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, numsweeps, numthreads);
else
  u = paric_async_mex(int32(Si)-1, int32(Sj)-1, Sa, u, numsweeps, numthreads);
end



