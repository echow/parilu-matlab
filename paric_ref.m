function [u resid] = paric_ref(a, u0, numsweeps)
% [u resid] = paric_ref(a, u0, numsweeps)
%   Synchronous updates
%
% a         = sparse matrix should be scaled to have unit diagonal
% u0        = input guess; pattern of u0 can be unrelated to matrix a
% numsweeps = number of nonlinear fixed-point sweeps
% u         = output computed upper triangular IC factor
% resid     = nonlinear residual norm history (optional)

if nargout > 1
    resid = zeros(numsweeps+1, 1);
    resid(1) = norm((a-u0'*u0).*spones(a),'fro');
end

% find nonzeros of u in column-major ordering
pat = spones(u0);
[Si Sj Sa] = find(a.*pat + eps*pat);
m = length(Sa);

% allocate output
u = u0;

for iter=1:numsweeps

    for k=1:m
        i = Si(k);
        j = Sj(k);
        s = Sa(k) - u0(1:i-1,i)'*u0(1:i-1,j);

        if (i ~= j)
            u(i,j) = s/u0(i,i);
        else
            if (s <= 0)
                fprintf('paric_ref: note pivot %f  sweep %d  index %d\n', s, iter, i);
            end
            u(i,j) = sqrt(s);
        end
    end

    u0 = u; % frozen u
    
    if nargout > 1
        resid(iter+1) = norm((a-u'*u).*spones(a),'fro');
    end

    % disp(nnz(imag(u)))
end
