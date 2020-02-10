function [l d u resid] = parildu_ref(a, l0, d0, u0, numsweeps)
% [l d u resid] = pariutdu_ref(a, l0, d0, u0, numsweeps)
% compute: l*d*u, where l and u have unit diagonal and d is diagonal
%   Synchronous updates
%
% a         = input sparse matrix
% l0        = input guess for L, with unit diagonal
% d0        = input guess for D, usually the diagonal of A
% u0        = input guess for U, with unit diagonal
% numsweeps = number of nonlinear fixed-point sweeps
% resid     = nonlinear residual norm history (optional)

if nargout > 3
    pat = spones(l0+u0);
    resid = zeros(numsweeps+1, 1);
    resid(1) = norm((a-l0*d0*u0).*spones(pat),'fro');
end

% find nonzeros of offdiagonal l (ordering doesn't matter in synchronous case)
pat = spones(l0) - speye(size(a));
[Si1 Sj1 Sa1] = find(a.*pat + eps*pat);
m1 = length(Sa1);

% find nonzeros of u (includes diagonal)
pat = spones(u0);
[Si2 Sj2 Sa2] = find(a.*pat + eps*pat);
m2 = length(Sa2);

% allocate space
l = l0;
u = u0;
dvec = diag(d0);
dvec0 = dvec;

for iter=1:numsweeps

    % lower triangular factor
    for k=1:m1
        i = Si1(k);
        j = Sj1(k);
        s = Sa1(k) - l0(i,1:j-1)*(dvec0(1:j-1).*u0(1:j-1,j));
        if (i == j), error('internal error'); end
        l(i,j) = s/dvec0(j);
    end

    % upper triangular factor and diagonal
    for k=1:m2
        i = Si2(k);
        j = Sj2(k);
        s = Sa2(k) - l0(i,1:i-1)*(dvec0(1:i-1).*u0(1:i-1,j));
        if (i ~= j)
            u(i,j) = s/dvec0(i);
        else
            dvec(i) = s;
        end
    end

    % reuse
    l0 = l;
    u0 = u;
    dvec0 = dvec;

    if nargout > 3
        resid(iter+1) = norm((a-l*diag(dvec)*u).*spones(pat),'fro');
    end
end
d = diag(dvec);

