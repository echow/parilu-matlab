function [u d resid] = parildlt_sync(a, u0, d0, numsweeps)
% [u d resid] = parildlt_sync(a, u0, numsweeps)
% compute: u'*d*u, where u has unit diagonal
% a = sparse matrix should be scaled to have unit diagonal
% u0 = input guess, where u should have unit diagonal
% d0 = input guess, usually the identity matrix
% numsweeps = number of nonlinear fixed-point sweeps
% resid = nonlinear residual norm history (Frobenius norm)
% note: more efficient to compute U because matrices are stored by columns

resid = zeros(numsweeps+1, 1);
resid(1) = norm((a-u0'*d0*u0).*spones(a),'fro');

% set of nonzeros of u in column-major ordering
pat = spones(u0);
[Si Sj Sa] = find(a.*pat + eps*pat);
m = length(Sa);

% allocate space
u = u0;
dvec = diag(d0);
dvec0 = dvec;

for iter=1:numsweeps

    for k=1:m
        i = Si(k);
        j = Sj(k);
        s = Sa(k) - u0(1:i-1,i)'*(dvec0(1:i-1).*u0(1:i-1,j));

        if (i ~= j)
            u(i,j) = s/dvec0(i);
        else
            if (s <= 0)
              fprintf('note pivot %f  sweep %d  index %d\n', s, iter, i);
            end
            dvec(i) = s;
        end
    end

    % reuse u0 and dvec0
    u0 = u;
    dvec0 = dvec;
    
    % nonlinear residual norm
    resid(iter+1) = norm((a-u'*diag(dvec)*u).*spones(a),'fro');
end
d = diag(dvec);
