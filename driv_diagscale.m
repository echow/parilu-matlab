function d = driv_diagscale

a60 = readmat(60);
a61 = readmat(61);

% smaller version of problem
p = 1:10002;
a60 = a60(p,p);
a61 = a61(p,p);

s60 = diagscale(a60);
s61 = diagscale(a61);

norm(a60-a61,'fro')/norm(a60,'fro')
norm(s60-s61,'fro')/norm(s60,'fro')

d60 = diag(a60);
d61 = diag(a61);

% if diagonal doesn't change, corresponding parts of matrix don't change
% that means scaling doesn't change...
d = d60-d61;
j = find(d==0);
nnz(j)
norm(a60(j,j)-a61(j,j),'fro')

spy(a60);
hold on
spy(a60-a61,'r');

function a = readmat(matnum)
filename = sprintf('/home/edmond/matrices.topopt/system_%05d.mat', matnum);
prob = load(filename);
a = prob.A;
a = (a+a')/2;

