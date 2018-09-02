function f=isOneInfinite(A)
n1=size(A,1);
n2=size(A,2);
f=false;
for i=1:n1
    for j=1:n2
        if (isinf(A(i,j)))
            f = true;
            return
        end
    end
end
end