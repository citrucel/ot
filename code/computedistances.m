function distanceMat = computedistances(allpoints)
n = size(allpoints,1);
distanceMat = zeros(n,n);
for i=1:n
    for j=i+1:n
        distanceMat(i,j) = pdist2(allpoints(i,:), allpoints(j,:),'euclidean');
        distanceMat(j,i) = distanceMat(i,j);
    end
end
end

