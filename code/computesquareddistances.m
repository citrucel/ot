function sqMinDistMat = computesquareddistances(allpoints, useMin, filterParameter)
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

