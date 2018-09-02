clc;
clear;

load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
global allpoints;
allpoints = data.allpoints;
dailyTemperatureStats = daytemperaturestats(allpoints);

global sqMinDistMat1;
sqMinDistMat1 = computesquareddistances(allpoints, true, 0);

global sqMinDistMat2;
sqMinDistMat2 = computesquareddistances(allpoints, true, 1);

global sqMinDistMat3;
sqMinDistMat3 = computesquareddistances(allpoints, true, 2);

tsneXY1 = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold1);
save('tsneXY1.mat', 'tsneXY1');
tsneXY2 = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold2);
save('tsneXY2.mat', 'tsneXY2');
tsneXY3 = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold3);
save('tsneXY3.mat', 'tsneXY3');

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
