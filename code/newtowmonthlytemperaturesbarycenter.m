clc;
clear;
load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
sampleSizes = data.sampleSizes;
allpoints = data.allpoints;
barycenter = barycenterndistributions(allpoints, sampleSizes)
sz = length(sampleSizes);
newSampleSizes = zeros(sz+1,1);
newSampleSizes(1:sz) = sampleSizes;
[mb,nb] = size(barycenter);
newSampleSizes(sz+1)=mb;
alldata.sampleSizes = newSampleSizes;
[m,n] = size(allpoints)
dayst = zeros(m + mb,n)
dayst(1:m,1:n) = allpoints(:,:)
dayst(m+1:m+mb,:)=barycenter(1:mb,:)
alldata.temperatures = dayst;
save('Adjusted2016Newtowntemperatures.mat','alldata')
