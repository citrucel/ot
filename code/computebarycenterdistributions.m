function baryCenterPoints=computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix)

baryCenter = 1:1:sampleSizes(1);
baryCenter = baryCenter';
count = 0;

for i=1:10
    newBaryCenter=computebarycenter(baryCenter, sampleSizes, distanceMatrix)
    tf = isequal(baryCenter,newBaryCenter);
    if tf == 1
        count = count + 1;
        if count == 1
            break
        end
    end
    baryCenter = newBaryCenter;
end
i
count
baryCenterPoints=findbarycenter(baryCenter, allpoints);
end


function extractedMatrix=extractdistancematrix(baryCenter, distanceMatrix)
n = size(baryCenter,1);
[~, nDistMatrix] = size(distanceMatrix);
extractedMatrix = zeros(n, nDistMatrix);
for i=1:n
    idx = baryCenter(i);
    extractedMatrix(i,:) = distanceMatrix(idx, :);
end
end

function baryCenter=computebarycenter(baryCenter,sampleSizes, distanceMatrix)

costMatrix = extractdistancematrix(baryCenter, distanceMatrix);

nbarycenter = size(baryCenter,1);
pbarycenter = ones(nbarycenter,1) * 1/nbarycenter;
nssz = size(sampleSizes,1);
n= 0;
nprev = 0;
cellMat=cell(1);
for i=1:nssz
    nsampleSizes = sampleSizes(i);
    n = n + nsampleSizes;
    p = ones(nsampleSizes,1) * 1/nsampleSizes;
    piMatrix = computewithlinprog(costMatrix(:, 1+nprev:n), pbarycenter, p);
    cellMat{1,i} = piMatrix;
    nprev = n;
end
aggPiMatrix=horzcat(cellMat{:});
alphaDistributionMat = aggPiMatrix * distanceMatrix;
minimalValues = min(alphaDistributionMat,[],2);
baryCenter=findalphaweightedclosestpoints(minimalValues, alphaDistributionMat);
end


function minIDX=findalphaweightedclosestpoints(minimalValues, alphaDistributionMat)
n = size(minimalValues,1);
minIDX = zeros(n,1);
for i=1:n
    minValue = minimalValues(i);
    [~,cols] = find(alphaDistributionMat==minValue);
    val = cols(1);
    minIDX(i)=val;
end
end

function baryCenter=findbarycenter(IDX, allpoints)
NIDX = unique(IDX);
%NIDX = IDX;
n = size(NIDX,1);
m = size(allpoints,2);
baryCenter = zeros(n,m);
for i=1:n
    idx = NIDX(i);
    baryCenter(i,1:m) = allpoints(idx,1:m);
end
end



