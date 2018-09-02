clc;
clear;
load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
sampleSizes = data.sampleSizes;
allpoints = data.allpoints;

distanceMatrix=computesquareddistances(allpoints, false, 0);
if isOneInfinite(distanceMatrix)
    return
end

nonMinBaryCenterPoints= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix);
nonMinDays=finddays(allpoints, nonMinBaryCenterPoints);

coldestDays = zeros(3,1);
warmestDays = zeros(3,1);
maxMinMonths = zeros(3,1);
meanMonths = zeros(3,1);
stdMonths = zeros(3,1);
hourlyDataSet = zeros(12,size(nonMinBaryCenterPoints,2));
[coldestDay, warmestDay,...
    maxMinMonth, meanMonth, stdMonth,...
    hourlyData] = computebarycenterstats(nonMinBaryCenterPoints);

coldestDays(1,1) = coldestDay;
warmestDays(1,1) = warmestDay;
maxMinMonths(1,1) = maxMinMonth;
meanMonths(1,1) = meanMonth;
stdMonths(1,1) = stdMonth;
hourlyDataSet(1:3,:) = hourlyData;

distanceMatrix2=computesquareddistances(allpoints, true, 0);
if isOneInfinite(distanceMatrix2)
    return
end
minBaryCenterPoints= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix2);
minDays=finddays(allpoints, minBaryCenterPoints);
[coldestDay, warmestDay,...
    maxMinMonth, meanMonth, stdMonth,...
    hourlyData] = computebarycenterstats(minBaryCenterPoints);
coldestDays(2,1) = coldestDay;
warmestDays(2,1) = warmestDay;
maxMinMonths(2,1) = maxMinMonth;
meanMonths(2,1) = meanMonth;
stdMonths(2,1) = stdMonth;
hourlyDataSet(4:6,:) = hourlyData;

distanceMatrix3=computesquareddistances(allpoints, true, 1);
if isOneInfinite(distanceMatrix3)
    return
end
minBaryCenterPoints2= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix3);
minDays2=finddays(allpoints, minBaryCenterPoints2);
[coldestDay, warmestDay,...
    maxMinMonth, meanMonth, stdMonth,...
    hourlyData] = computebarycenterstats(minBaryCenterPoints2);
coldestDays(3,1) = coldestDay;
warmestDays(3,1) = warmestDay;
maxMinMonths(3,1) = maxMinMonth;
meanMonths(3,1) = meanMonth;
stdMonths(3,1) = stdMonth;
hourlyDataSet(7:9,:) = hourlyData;

distanceMatrix4=computesquareddistances(allpoints, true, 2);
if isOneInfinite(distanceMatrix4)
    return
end
minBaryCenterPoints3= computebarycenterdistributions(sampleSizes, allpoints, distanceMatrix4);
minDays3=finddays(allpoints, minBaryCenterPoints3);
[coldestDay, warmestDay,...
    maxMinMonth, meanMonth, stdMonth,...
    hourlyData] = computebarycenterstats(minBaryCenterPoints3);
coldestDays(4,1) = coldestDay;
warmestDays(4,1) = warmestDay;
maxMinMonths(4,1) = maxMinMonth;
meanMonths(4,1) = meanMonth;
stdMonths(4,1) = stdMonth;
hourlyDataSet(10:12,:) = hourlyData;

hourlyData = computebyhour(allpoints);
hourlyDataSet(13:15,:) = hourlyData;

disp('Coldest Day on average orig-min')
coldestDays
disp('Warmest Day on average orig-min')
warmestDays
disp('Max-Min Temperature on average orig-min')
maxMinMonths
disp('Mean orig-min')
meanMonths
disp('Std orig-min')
stdMonths
plothourlytemperatures(hourlyDataSet);

function plothourlytemperatures(hourlyData)
line(1:24, hourlyData(2,:), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);
% Plot upper and lower bounds
l1=line(1:24, hourlyData(1,:), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);
l2=line(1:24, hourlyData(3,:), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 3);

line(1:24, hourlyData(5,:), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 3);
% Plot upper and lower bounds
l3=line(1:24, hourlyData(4,:), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 3);
l4=line(1:24, hourlyData(6,:), 'Color', 'g', 'LineStyle', '--', 'LineWidth', 3);

line(1:24, hourlyData(8,:), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 3);
% Plot upper and lower bounds
l5=line(1:24, hourlyData(7,:), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 3);
l6=line(1:24, hourlyData(9,:), 'Color', 'b', 'LineStyle', '--', 'LineWidth', 3);


line(1:24, hourlyData(11,:), 'Color', 'm', 'LineStyle', '-', 'LineWidth', 3);
% Plot upper and lower bounds
l7=line(1:24, hourlyData(12,:), 'Color', 'm', 'LineStyle', '-', 'LineWidth', 3);
l8=line(1:24, hourlyData(10,:), 'Color', 'm', 'LineStyle', '-', 'LineWidth', 3);

line(1:24, hourlyData(14,:), 'Color', 'c', 'LineStyle', '-', 'LineWidth', 3);
% Plot upper and lower bounds
l9=line(1:24, hourlyData(15,:), 'Color', 'c', 'LineStyle', '-', 'LineWidth', 3);
l10=line(1:24, hourlyData(13,:), 'Color', 'c', 'LineStyle', '-', 'LineWidth', 3);

set(get(get(l1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l8,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l9,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(l10,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xtickat = 0:1:length(hourlyData(2,:))-1;  %every hour from 00:00
set(gca, 'XTick', xtickat, 'XTickLabel', cellstr( num2str( mod(round(xtickat .' ./ 1),24) ) ) )

xlabel('Hour of the Day')
ylabel('Hourly Temperature')
axis tight;
legend('Location','southeast')
legend('Init.Matrix','Min.Matrix.Threshold 1','Min.Matrix.Threshold 2','Min.Matrix.d-distance','Init.Data')

end

function days=finddays(allpoints, barycenter)
[n2,~] = size(barycenter);
days = zeros(n2,1);
for i=1:n2
    [~, index]=ismember(barycenter(i,:),allpoints,'rows');
    days(i)=index;
end
end


