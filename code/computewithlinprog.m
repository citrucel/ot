function pimatrix=computewithlinprog(F, p0, p1)
n0 = length(p0);
n1 = length(p1);
flat = @(x)x(:);
Cols = @(n0,n1)sparse( flat(repmat(1:n1, [n0 1])), ...
    flat(reshape(1:n0*n1,n0,n1) ), ...
    ones(n0*n1,1) );
Rows = @(n0,n1)sparse( flat(repmat(1:n0, [n1 1])), ...
    flat(reshape(1:n0*n1,n0,n1)' ), ...
    ones(n0*n1,1) );
Sigma = @(n0,n1)[Rows(n0,n1);Cols(n0,n1)];
lb = zeros(n0*n1, 1);
ub = ones(n0*n1, 1);
%linprog(f,A,b,Aeq,beq,lb,ub)
otransp = @(F,p0,p1)reshape( linprog(F(:),[],[],...
    Sigma(n0,n1), [p0(:);p1(:)], lb, ub), [n0 n1]);
pimatrix = otransp(F,p0,p1);
end

