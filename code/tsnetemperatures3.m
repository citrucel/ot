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
%'MarkerFaceColor',c(j,:),'Marker',m{i})
scatterPlot('Distance Threshold 1', tsneXY1(:,1),tsneXY1(:,2), dailyTemperatureStats, idxSeasons)

h(5)=subplot(3,2,5);
scatterPlot('Distance Threshold 2', tsneXY2(:,1),tsneXY2(:,2), dailyTemperatureStats, idxSeasons)

h(6)=subplot(3,2,6);
scatterPlot('d-distance', tsneXY3(:,1),tsneXY3(:,2), dailyTemperatureStats, idxSeasons)

pos = get(h,'position');
set(h(1),'position',pos{1});
set(h(2),'position',pos{3});
set(h(3),'position',pos{5});
set(h(4),'position',pos{2});
set(h(5),'position',pos{4});
set(h(6),'position',pos{6});

colorbar

function scatterPlot(titleStr, x, y, dailyTemperatureStats, seasons)
[winter,winterTemperatureStats,spring,springTemperatureStats, ...
    summer,summerTemperatureStats, fall, fallTemperatureStats]=partitionBySeason(x, y, seasons, dailyTemperatureStats);
scatter(winter(:,1),winter(:,2),winterTemperatureStats(:,1),winterTemperatureStats(:,2), 'filled', 'p')
hold on
scatter(spring(:,1),spring(:,2),springTemperatureStats(:,1),springTemperatureStats(:,2), 'filled', 's')
hold on
scatter(summer(:,1),summer(:,2),summerTemperatureStats(:,1),summerTemperatureStats(:,2), 'filled', 'o')
hold on
scatter(fall(:,1),fall(:,2),fallTemperatureStats(:,1),fallTemperatureStats(:,2), 'filled', 'd')
hold on
lgd=legend('winter','spring','summer','fall','Location','northeast');
lgd.FontSize = 10;
lgd.TextColor = 'black';
title(titleStr)
end

function [winter,winterTemperatureStats,spring,springTemperatureStats, summer,summerTemperatureStats, ...
    fall, fallTemperatureStats]=... 
    partitionBySeason(x, y, seasons, dailyTemperatureStats)
n = size(x,1);
% winter = zeros(n,2);
% winterTemperatureStats = zeros(n,2);
% spring = zeros(n,2);
% springTemperatureStats = zeros(n,2);
% summer = zeros(n,2);
% summerTemperatureStats = zeros(n,2);
% fall = zeros(n,2);
% fallTemperatureStats = zeros(n,2);
j1 = 1;
j2 = 1;
j3 = 1;
j4 = 1;
K = 3;
for i=1:n
    % Winter
    if seasons(i) == 1
        winter(j1,1) = x(i);
        winter(j1,2) = y(i);
        winterTemperatureStats(j1,1) = K * dailyTemperatureStats(i,1);
        winterTemperatureStats(j1,2) = dailyTemperatureStats(i,2);
        j1 = j1 +1;
        % Spring
    elseif seasons(i) == 2
        spring(j2,1) = x(i);
        spring(j2,2) = y(i);
        springTemperatureStats(j2,1) = K * dailyTemperatureStats(i,1);
        springTemperatureStats(j2,2) = dailyTemperatureStats(i,2);
        j2 = j2 +1;
        % Summer
    elseif seasons(i) == 3
        summer(j3,1) = x(i);
        summer(j3,2) = y(i);
        summerTemperatureStats(j3,1) = K * dailyTemperatureStats(i,1);
        summerTemperatureStats(j3,2) = dailyTemperatureStats(i,2);
        j3 = j3 +1;
        % Fall
    elseif seasons(i) == 4
        fall(j4,1) = x(i);
        fall(j4,2) = y(i);
        fallTemperatureStats(j4,1) = K * dailyTemperatureStats(i,1);
        fallTemperatureStats(j4,2) = dailyTemperatureStats(i,2);
        j4 = j4 +1;
    end
end

% winter = filterZeros(winter);
% winterTemperatureStats = filterZeros(winterTemperatureStats);
% spring = filterZeros(spring);
% springTemperatureStats = filterZeros(springTemperatureStats);
% summer = filterZeros(summer);
% summerTemperatureStats = filterZeros(summerTemperatureStats);
% fall = filterZeros(fall);
% fallTemperatureStats = filterZeros(fallTemperatureStats);

end


function nv = filterZeros(v)
x  = v(:,1);
nvx = x(x ~= 0);
y  = v(:,2);
nvy = y(y ~= 0);
nv = zeros(length(nvx),2);
nv(:,1) =nvx;
nv(:,2) =nvy;
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



