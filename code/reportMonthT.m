clc;
clear;

load('Adjusted2016Newtowntemperatures.mat', 'alldata');
sampleSizes = alldata.sampleSizes;
[ms,ns] = size(sampleSizes);
temperatures = alldata.temperatures;
[m,n] = size(temperatures);
[meanmonths, stdmonths] = computemeanstd(sampleSizes, temperatures);
[coldestMonths, warmestMonths, maxMinMonths] = computemonthlystats(sampleSizes, temperatures);
hourlyData = computebyhour(temperatures(1:m-sampleSizes(ms),:));
barycenter = temperatures(m-sampleSizes(ms)+1:m, :);
hourlyDataBarycenter = computebyhour(barycenter);
[minT, maxT] = computeminmaxhourly(sampleSizes, temperatures);

plotData(meanmonths, stdmonths,...
    coldestMonths, warmestMonths, maxMinMonths,...
    hourlyData, hourlyDataBarycenter, minT, maxT)

function [meanMonths, stdMonths] = computemeanstd(sampleSizes, temperatures)
[m,~] = size(sampleSizes);
meanMonths = zeros(m,1);
stdMonths = zeros(m,1);
sz = 0;
for i=1:m
    nsz = sz + sampleSizes(i);
    temperatureMonths = temperatures(sz+1:nsz, :);
    temperatureMonths = temperatureMonths(:);
    meanMonths(i)= mean(temperatureMonths);
    stdMonths(i)=std(temperatureMonths,1);
    sz = nsz;
end
end

function [coldestMonths, warmestMonths, maxMinMonths] = computemonthlystats(sampleSizes, temperatures)
[m,~] = size(sampleSizes);
coldestMonths = zeros(m,1);
warmestMonths = zeros(m,1);
maxMinMonths = zeros(m,1);
sz = 0;
for i=1:m
    nsz = sz + sampleSizes(i);
    temperatureMonths = temperatures(sz+1:nsz, :);
    [coldestDay,warmestDay] = coldestwarmestdays(temperatureMonths);
    coldestMonths(i) = coldestDay;
    warmestMonths(i) = warmestDay;
    maxMinT = variabilitydays(temperatureMonths);
    maxMinMonths(i) = maxMinT;
    sz = nsz;
end
end

function plotData(meanmonths, stdmonths,...
    coldestMonths, warmestMonths, maxMinMonths,...
    hourlyData, hourlyDataBarycenter,...
    minT, maxT)

Figure1=figure;
set(Figure1,'defaulttextinterpreter','latex');
h(1) = subplot(3,2,1);
plotmonthlytemperatures(h(1), coldestMonths, 'Coldest T')
h(2) = subplot(3,2,2);
plotmonthlytemperatures(h(2), warmestMonths, 'Warmest T')
h(3) = subplot(3,2,3);
plotmonthlytemperatures(h(3), maxMinMonths, 'Variability')
h(4) = subplot(3,2,4);
plotmonthlymeantemperatures(h(4), meanmonths, stdmonths)
h(5) = subplot(3,2,5);
plothourlytemperatures(hourlyData, hourlyDataBarycenter)
h(6) = subplot(3,2,6);
plotminTmaxT(h(6), minT, maxT)
end

function plothourlytemperatures(hourlyData, hourlyDataBarycenter)
l1=line(1:24, hourlyData(2,:), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
% Plot upper and lower bounds
l2=line(1:24, hourlyData(1,:), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
l3=line(1:24, hourlyData(3,:), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

line(1:24, hourlyDataBarycenter(2,:), 'Color', 'g', 'LineStyle', '-', 'LineWidth', 2);
% Plot upper and lower bounds
l4=line(1:24, hourlyDataBarycenter(1,:), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
l5=line(1:24, hourlyDataBarycenter(3,:), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);

set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xtickat = 0:1:length(hourlyData(2,:))-1;  %every hour from 00:00
set(gca, 'XTick', xtickat, 'XTickLabel', cellstr( num2str( mod(round(xtickat .' ./ 1),24) ) ) )

xlabel('Hour of the Day')
ylabel('Hourly Temperature')
axis tight;

end

function plotmonthlymeantemperatures(h, data, stddata)
l1=line(1:12, data(1:12), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
datastd = data(1:12) - stddata(1:12);
% Plot upper and lower bounds
l2=line(1:12, datastd, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
datastd = data(1:12) + stddata(1:12);
l3=line(1:12, datastd, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

repdata = repmat(data(13),12,1);
line(1:12, repmat(data(13),12,1), 'Color', 'g', 'LineStyle', '-', 'LineWidth', 2);
% Plot upper and lower bounds
repstddata = repmat(stddata(13),12,1);
datastd = repdata - repstddata;
l4=line(1:12, datastd, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
datastd = repdata + repstddata;
l5=line(1:12, datastd, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);

set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


set(h,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})

xlabel('2016')
ylabel('Mean +/- $\sigma$ ')
axis tight

end

function plotmonthlytemperatures(h, data, ylabelText)
l1=line(1:12, data(1:12), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
line(1:12, repmat(data(13),12,1), 'Color', 'g','LineStyle','-','LineWidth',2)
set(h,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('2016')
ylabel(ylabelText)
axis tight
end

function plotminTmaxT(h, minT, maxT)
j = size(minT,2);
for i=1:j
    plotminmax(minT(:,i), maxT(:,i));
end
set(h,'xtick',1:12,...
    'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
xlabel('2016')
ylabel('T [4am,8am,12pm,16pm]','fontweight','bold')
axis tight
end

function plotminmax(min, max)
line(1:12, min(1:12), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
line(1:12, repmat(min(13),12,1),'Color', 'g','LineStyle','--','LineWidth',2)
line(1:12, max(1:12), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
line(1:12, repmat(max(13),12,1),'Color', 'g','LineStyle','-','LineWidth',2)
end

