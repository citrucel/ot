clc;
clear;
load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
sampleSizes = data.sampleSizes;
allpoints = data.allpoints;
distanceMatrix=computedistances(allpoints, true, 2);
if isOneInfinite(distanceMatrix)
    return
end
minBaryCenterPoints= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix);

function sqMinDistMat = computedistances(allpoints, useMin, filterParameter)
n = size(allpoints,1);
distanceMat = zeros(n,n);
for i=1:n
    for j=i+1:n
        distanceMat(i,j) =pdist2(allpoints(i,:), allpoints(j,:),'euclidean');
        distanceMat(j,i) = distanceMat(i,j);
    end
end
if useMin == true
    disp('Minimal distance matrix')
    minDistMat= mindist(distanceMat, filterParameter);
    disp('Squared minimal distance matrix')
    sqMinDistMat = (minDistMat.^2)./2;
else
    disp('Squared minimal distance matrix')
    sqMinDistMat = (distanceMat.^2)./2;
end
end

function f=isOneInfinite(A)
[n1,n2]=size(A);
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