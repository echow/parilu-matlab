function [l u] = parilu(a, l, u, numsweeps, numthreads)
% [l u] = parilu(a, l, u, numsweeps, numthreads)
%  Input l and u are initial guesses.
%  Sparsity pattern of l and u are preserved on output.
%  
% If you want l and u to contain structural zero elements,
% then you must put in small values.

% following code is very slow
%
% [Si Sj] = find(spones(l)+spones(u));
% m = length(Si);
% Sa = zeros(m,1);
% for k=1:m
%   Sa(k) = a(Si(k),Sj(k));
% end

pat = spones(l-diag(diag(l))+u);
[Si Sj Sa] = find(a.*pat + 1e-99*pat);

% convert indices Si and Sj to 0-based
l = l';
[l u] = parilu_mex(int32(Si)-1, int32(Sj)-1, Sa, l, u, numsweeps, numthreads);
l = l';


