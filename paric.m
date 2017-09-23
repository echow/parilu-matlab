function u = paric(a, u, numsweeps, numthreads)
% u = paric(a, u, numsweeps, numthreads)

% notes:
% ordering for find is important
% 0-based indexing for mex file

pat = spones(u);
[Si Sj Sa] = find(a.*pat + eps*pat);
%u = paric_sync_mex(int32(Si)-1, int32(Sj)-1, Sa, u, numsweeps, numthreads);
 u = paric_async_mex(int32(Si)-1, int32(Sj)-1, Sa, u, numsweeps, numthreads);



