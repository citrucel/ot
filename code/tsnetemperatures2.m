clc;
clear;

load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
allpoints = data.allpoints;

load('Adjusted2016Newtowntemperatures.mat', 'alldata');
sampleSizes = alldata.sampleSizes;
[ms,ns] = size(sampleSizes);
temperatures = alldata.temperatures;
[m,n] = size(temperatures);
coldwarmlabels=coldwarmlabel(allpoints , sampleSizes, temperatures);
coldwarmCategories=categorical(coldwarmlabels,{'Cold','Mild','Hot '},'Ordinal',true);

tsneXY1 = load('tsneXY1.mat','tsneXY1');
tsneXY1 = tsneXY1.tsneXY1;
tsneXY2 = load('tsneXY2.mat','tsneXY2');
tsneXY2 = tsneXY2.tsneXY2;
tsneXY3 = load('tsneXY3.mat','tsneXY3');
tsneXY3 = tsneXY3.tsneXY3;

h(1)=subplot(3,2,1);
clrs = ['bgy'];
gscatter(tsneXY1(:,1),tsneXY1(:,2),coldwarmCategories,clrs)
title('Distance Threshold 1')

h(2)=subplot(3,2,2);
gscatter(tsneXY2(:,1),tsneXY2(:,2),coldwarmCategories,clrs)
title('Distance Threshold 2')

h(3)=subplot(3,2,3);
gscatter(tsneXY3(:,1),tsneXY3(:,2),coldwarmCategories,clrs)
title('d-distance')

seasons = data.season';
idxSeasons = double(seasons) ;
dailyTemperatureStats = daytemperaturestats(allpoints);

h(4)=subplot(3,2,4);
scatter(tsneXY1(:,1),tsneXY1(:,2), dailyTemperatureStats(:,1),dailyTemperatureStats(:,2))
title('Distance Threshold 1')

h(5)=subplot(3,2,5);
scatter(tsneXY2(:,1),tsneXY2(:,2), dailyTemperatureStats(:,1),dailyTemperatureStats(:,2))
title('Distance Threshold 2')

h(6)=subplot(3,2,6);
scatter(tsneXY3(:,1),tsneXY3(:,2), dailyTemperatureStats(:,1),dailyTemperatureStats(:,2))
title('d-distance')

pos = get(h,'position');
set(h(1),'position',pos{1});
set(h(2),'position',pos{3});
set(h(3),'position',pos{5});
set(h(4),'position',pos{2});
set(h(5),'position',pos{4});
set(h(6),'position',pos{6});

colorbar

function dailyTemperatureStats = daytemperaturestats(allpoints)
m = size(allpoints,1);
dailyTemperatureStats = zeros(m,2);
for i=1:m
    dailyTemperatureStats(i,1) = max(allpoints(i,:)) - min(allpoints(i,:)); 
    dailyTemperatureStats(i,2) = mean(allpoints(i,:));
end
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



