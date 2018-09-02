clc;
clear;

load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
global allpoints;
allpoints = data.allpoints;
seasons = data.season';

global sqMinDistMat1;
sqMinDistMat1 = computesquareddistances(allpoints, true, 0);

%Y = tsne(allpoints,'Algorithm','barneshut','Distance',@distancethrehsold1);
Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold1);
h(1)=subplot(3,2,1);
gscatter(Y(:,1),Y(:,2),seasons)
title('Distance Threshold1')

global sqMinDistMat2;
sqMinDistMat2 = computesquareddistances(allpoints, true, 1);
Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold2);
h(2)=subplot(3,2,2);
gscatter(Y(:,1),Y(:,2),seasons)
title('Distance Threshold2')

global sqMinDistMat3;
sqMinDistMat3 = computesquareddistances(allpoints, true, 2);
Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold3);
h(3)=subplot(3,2,3);
gscatter(Y(:,1),Y(:,2),seasons)
title('d-distance')

load('Adjusted2016Newtowntemperatures.mat', 'alldata');
sampleSizes = alldata.sampleSizes;
[ms,ns] = size(sampleSizes);
temperatures = alldata.temperatures;
[m,n] = size(temperatures);
coldwarmlabels=coldwarmlabel(allpoints , sampleSizes, temperatures);
coldwarmCategories=categorical(coldwarmlabels,{'Cold','Mild','Hot '},'Ordinal',true);

Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold1);
h(4)=subplot(3,2,4);
gscatter(Y(:,1),Y(:,2),coldwarmCategories)
title('Distance Threshold1')

Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold2);
h(5)=subplot(3,2,5);
gscatter(Y(:,1),Y(:,2),coldwarmCategories)
title('Distance Threshold2')


Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold3);
h(6)=subplot(3,2,6);
gscatter(Y(:,1),Y(:,2),coldwarmCategories)
title('d-distance')

pos = get(h,'position');
set(h(1),'position',pos{1});
set(h(2),'position',pos{3});
set(h(3),'position',pos{5});
set(h(4),'position',pos{2});
set(h(5),'position',pos{4});
set(h(6),'position',pos{6});

function d = distancethrehsold1(x, y)
global allpoints sqMinDistMat1
i =find(ismember(allpoints,x,'rows'));
j =find(ismember(allpoints,y,'rows'));
d = sqMinDistMat1(i,j)';
end

function d = distancethrehsold2(x, y)
global allpoints sqMinDistMat2
i =find(ismember(allpoints,x,'rows'));
j =find(ismember(allpoints,y,'rows'));
d = sqMinDistMat2(i,j)';
end

function d = distancethrehsold3(x, y)
global allpoints sqMinDistMat3
i =find(ismember(allpoints,x,'rows'));
j =find(ismember(allpoints,y,'rows'));
d = sqMinDistMat3(i,j)';
end

function labels=coldwarmlabel(allpoints , sampleSizes, temperatures)
[coldestMonths, warmestMonths] = computemonthlystats(sampleSizes, temperatures);
m = size(allpoints,1);
labels{m} = [];
ms = length(sampleSizes);
coldBarycenter = coldestMonths(ms);
warmBarycenter = warmestMonths(ms);
for i=1:m
    avg = mean(allpoints(i,:));
    if (avg < coldBarycenter)
        labels{i} = 'Cold'; 
    elseif (avg > warmBarycenter)
        labels{i} = 'Hot ';
    else
        labels{i} = 'Mild';        
    end
end
end

function [coldestMonths, warmestMonths] = computemonthlystats(sampleSizes, temperatures)
[m,~] = size(sampleSizes);
coldestMonths = zeros(m,1);
warmestMonths = zeros(m,1);
sz = 0;
for i=1:m
    nsz = sz + sampleSizes(i);
    temperatureMonths = temperatures(sz+1:nsz, :);
    [coldestDay,warmestDay] = coldestwarmestdays(temperatureMonths);
    coldestMonths(i) = coldestDay;
    warmestMonths(i) = warmestDay;
    sz = nsz;
end
end
