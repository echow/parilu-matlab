function [l u resid] = parilu_ref(a, l, u, numsweeps)
% [l u resid] = parilu_ref(a, l, u, numsweeps)
%   Synchronous updates
%
% a         = input sparse matrix
% l, u      = input guesses; patterns can be unrelated to matrix a;
%             diagonal of l should be all ones; output factors
% numsweeps = number of nonlinear fixed-point sweeps
% resid     = nonlinear residual norm history (optional)

if nargout > 2
    pat = spones(spones(l)+spones(u));
    resid = zeros(numsweeps+1, 1);
    resid(1) = norm((a-l*u).*spones(pat),'fro');
end

% find nonzeros of u in column-major ordering
pat = spones(spones(l) + spones(u));
[Si Sj Sa] = find(a.*pat + eps*pat);
m = length(Sa);

for iter=1:numsweeps

    l0 = l'; % frozen l in transposed form
    u0 = u;  % frozen u
    
    for k=1:m
        i = Si(k);
        j = Sj(k);
        if (i>j)
            l(i,j) = ( Sa(k) - l0(:,i)'*u0(:,j) + l0(j,i)*u0(j,j) ) / u0(j,j);
        else
            u(i,j) =   Sa(k) - l0(:,i)'*u0(:,j) + l0(i,i)*u0(i,j);
            % note l(i,i) is 1
        end
    end
    
    if nargout > 2
        resid(iter+1) = norm((a-l*u).*spones(pat),'fro');
    end
end

% below is the version if l is not transposed
%
% for iter=1:numsweeps
% 
%     l0 = l; % frozen l
%     u0 = u; % frozen u
% 
%     for k=1:m
%         i = Si(k);
%         j = Sj(k);
%         if (i>j)
%             l(i,j) = Sa(k) - l0(i,:)*u0(:,j) + l0(i,j)*u0(j,j);
%             l(i,j) = l(i,j) / u0(j,j);
%         else
%             u(i,j) = Sa(k) - l0(i,:)*u0(:,j) + l0(i,i)*u0(i,j);
%         end
%     end
% 
%     resid(iter+1) = norm((a-l*u).*spones(a),'fro');
% end

