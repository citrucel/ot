clc;
clear;


n0 = 500
n1 = 500
n = n0 + n1;
xy = zeros(n,2);
xyd = createnoisydistribution(n0, 0, 0.8);
xy(1:n0,1) = xyd(:,1);
xy(1:n0,2) = xyd(:,2);
disp('Mu1 distribution')
xy(1:n0,:)
xyd = createnoisydistribution(n1, pi/2, 0.8);
xy(n0+1:n,1) = xyd(:,1);
xy(n0+1:n,2) = xyd(:,2);
disp('Mu2 distribution')
xy(n0+1:n,:)

xymu1 = computeexpdistribution(n0, n1, xy, true, 0);
xymu2 = computeexpdistribution(n0, n1, xy, true, 1);
xymu3 = computeexpdistribution(n0, n1, xy, true, 2);
plotresults(n0, xy, xymu1, xymu2, xymu3)

function xymu = computeexpdistribution(n0, n1, xy, useMin, filterParameter)
sqMinDistMat = computesquareddistances(xy, useMin, filterParameter);
n = n0 + n1;
minCost = sqMinDistMat(1:n0, n0+1:n);
if isOneInfinite(minCost)
    return
end

p0 = ones(n0,1) * 1/n0;
p1 = ones(n1,1) * 1/n1;
pimatrix = computewithlinprog(minCost, p0, p1);
xymu = (n0 .* pimatrix) * xy(n0+1:n, :);
end

function xy=createdistribution(n, theta0, sigma)
theta = theta0 + sigma.*randn(1,n);
xy = zeros(n,2);
xy(1:n,1) = cos(theta);
xy(1:n,2) = sin(theta);
end


function xy=createnoisydistribution(n, theta0, sigma)
theta = theta0 + sigma.*randn(1,n);
xy = zeros(n,2);
noise = rand( 1, n) * .1;
r = 1 + noise;
xy(1:n,1) = r .* cos(theta);
xy(1:n,2) = r .* sin(theta);
end

function plotresults(n0, xy, xymu1, xymu2, xymu3)
n = size(xy,1);
pl1 = plotpoints(xy(1:n0,1),xy(1:n0,2),'r');
set(get(get(pl1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
pl2 = plotpoints(xymu1(:,1), xymu1(:,2),'g');
set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
pl3 = plotpoints(xymu2(:,1), xymu2(:,2),'b');
set(get(get(pl3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
pl4 = plotpoints(xymu3(:,1), xymu3(:,2),'m');
set(get(get(pl4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold on
for k=1:n0
    plot([xy(k,1) xymu1(k,1)], [xy(k,2) xymu1(k,2)], 'g' , 'LineStyle', '--', 'LineWidth', 0.5);
    plot([xy(k,1) xymu2(k,1)], [xy(k,2) xymu2(k,2)], 'b' , 'LineStyle', '--', 'LineWidth', 0.5);
    plot([xy(k,1) xymu3(k,1)], [xy(k,2) xymu3(k,2)], 'm' , 'LineStyle', '--', 'LineWidth', 0.5);
end
axis square
lgd=legend('Min.Matrix.1-distance','Min.Matrix.2-distance','Min.Matrix.d-distance','Location','southwest');
lgd.FontSize = 10;
lgd.TextColor = 'black';
% txt = texlabel(['|X|=' num2str(n0) ' and |Y|=' num2str(n-n0)]);
% title(txt, 'FontSize', 14)
end


function pl = plotpoints(x, y, col)
myplot = @(x,y,col)plot(x,y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);
pl = myplot(x, y, col);
end


