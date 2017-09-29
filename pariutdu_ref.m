function [u d resid] = pariutdu_ref(a, u0, d0, numsweeps)
% [u d resid] = pariutdu_ref(a, u0, d0, numsweeps)
% compute: u'*d*u, where u has unit diagonal and d is diagonal
%   Synchronous updates
%
% a         = sparse matrix should be scaled to have unit diagonal
% u0        = input guess, with unit diagonal
% d0        = input guess, usually the identity matrix, but may contain
%             both positive and negative entries
% numsweeps = number of nonlinear fixed-point sweeps
% resid     = nonlinear residual norm history (optional)

% note: more efficient to compute U because matrices are stored by columns

if nargout > 2
    resid = zeros(numsweeps+1, 1);
    resid(1) = norm((a-u0'*d0*u0).*spones(a),'fro');
end

% find nonzeros of u in column-major ordering
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
    
    if nargout > 2
        resid(iter+1) = norm((a-u'*diag(dvec)*u).*spones(a),'fro');
    end
end
d = diag(dvec);
