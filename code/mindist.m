function D=mindist(A, filterParameter)
% A: data
% filterParameter: how to filter the initial matrix

if filterParameter == 0
    R=filtermatrixabovethreshold(A);
elseif filterParameter == 1
    R=filtermatrix(A);    
elseif filterParameter == 2
    R = power(A);
else
    R = A;
end

D=computemindistuntilnoInf(A, R);
end


function DM=computemindistuntilnoInf(A, D) 
    loopcount = 20;
    hasInf = true;
    while loopcount ~= 0 && hasInf
        DM = computemindist(D);
        hasInf = isOneInfinite(DM);
        if hasInf
            loopcount = loopcount - 1;
            D=replaceInf(A, DM);
        end
    end
end


% min_dist: Floyd-Warshall algorithm
function B = computemindist(A)
	n = size(A, 1);
    B = A(:,:);
    for k=1:n
        for i=1:n
            for j=1:n
                s = B(i,k) + B(k,j);
                if B(i,j) > s
                    %fprintf('i=%d,j=%d,k=%d, iv=%d, v=%d\n', i,j,k, B(i,j) ,s)
                    B(i,j) = s;
                end
            end
        end
    end
end

% function threshold = maxmin(A)
% 	n = size(A, 1);
%     minv = inf(n,1);
%     for i=1:n
%         minval = inf;
%         for j=1:n
%             if i ~= j
%                 if A(i,j) < minval
%                     minval = A(i,j);
%                 end
%             end
%         end
%         minv(i)=minval;
%     end
%     threshold = max(minv);
% end

function new=filtermatrixabovethreshold(A)
    threshold = maxmin(A);
    %threshold = min(weigth .* median(A));
    %threshold = 2*p/np;
    fprintf('threshold=%g\n', threshold)
    new = A;
    new(abs(new)>threshold) = inf;
end

function threshold = maxmin(A)
    idx = 1==eye(size(A));
    new = A;
    new(idx) = inf;
    threshold =  max(min(new));
end


function new=filtermatrix(A)
    [n,m]=size(A);
    new = inf(n, m);
    idx = 1==eye(size(A));
    new(idx) = 0;
    B = A;
    B(idx) = inf;
    [minv,idx] = min(B, [], 1);
    for i=1:n
        j = idx(i);
        new(i, j) = minv(i);
        new(j,i) = new(i,j);
    end
end

function new=power(A)
    [n,m]=size(A);
    new = A^2;
    idx = 1==eye(size(A));
    new(idx) = 0;
end

function C=replaceInf(A, B)
    n=size(B,1);
    C=zeros(n,n);
    for i=1:n
      [mv, jmv] = minv(i, A(i,:), B(i,:));
      if (isOneInfinite(B(i,:)))
         C(i, jmv) = mv;
         C(jmv,i) = C(i,jmv);
      end
      for j=1:n
          if j ~= jmv 
              C(i,j)= B(i,j);
              C(j,i) = C(i,j);
          end
      end
    end
end


function [minval, minj]=minv(i, v1, v2)
  minval = Inf;
  minj = i;
  n=length(v1);
  for j=1:n
    if i ~= j
        if isinf(v2(j)) && v1(j) < minval
            minval = v1(j);
            minj = j;
        end
    end
  end
end


function checkResult(A,B,T,threshold)
    n = size(A, 1);
    for i=1:n
        for j=1:n
            if T(i,j) == 0 && B(i,j) > threshold && B(i,j) == A(i,j)
                fprintf('i=%d,j=%d,B(i,j)=%g\n', i,j,B(i,j))
            end
        end
    end
end