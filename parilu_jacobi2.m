function [l u resid] = parilu_jacobi(a, l, u, numsweeps)
% a = sparse matrix should be scaled to have unit diagonal
% l u = input guess; output computed ILU(0) factors with unit triangular l
% numsweeps = number of nonlinear fixed-point sweeps
% resid = residual norm history (Frobenius norm)

resid = zeros(numsweeps+1, 1);
resid(1) = norm((a-l*u).*spones(a),'fro');

% set of nonzeros of a in column-major ordering
[Si Sj] = find(a);
m = length(Si);
Sa = zeros(m,1);
for k=1:m
  Sa(k) = a(Si(k),Sj(k));
end

for iter=1:numsweeps

    l0 = l; % frozen l
    u0 = u; % frozen u
    
    for k=1:m
        i = Si(k);
        j = Sj(k);
        if (i>j)
            l(i,j) = Sa(k) - l0(i,:)*u0(:,j) + l0(i,j)*u0(j,j);
            l(i,j) = l(i,j) / u0(j,j);
        else
            u(i,j) = Sa(k) - l0(i,:)*u0(:,j) + l0(i,i)*u0(i,j);
            % note note l(i,i) is 1
        end
    end
    
    % compute residual norm
    xxx = (a-l*u).*spones(a);
    imagesc(log(abs(xxx)));
    colorbar
    pause
    resid(iter+1) = norm((a-l*u).*spones(a),'fro');
end
