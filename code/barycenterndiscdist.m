
clear
clc

n1 = 60;
disp('Mu1 distribution')
mu1 = createdistribution(n1,0,0.2)

n2 = 80;
disp('Mu2 distribution')
mu2 = createdistribution(n2,pi/2,0.2)
allpoints = vertcat(mu1, mu2)
sampleSizes = [n1 n2]';
baryCenterPoints=barycenterndistributions(allpoints, sampleSizes)

sampleShapes = ['s';'d'];
sampleColors = ['y';'g'];
h = figure;
h1 = subplot(1,3,1);
titletext= "Two distributions";
plotresults(sampleSizes, sampleShapes, sampleColors,...
    allpoints, baryCenterPoints, titletext)

n1 = 60;
disp('Mu1 distribution')
mu1 = createdistribution(n1,0,0.2)
n2 = 80;
disp('Mu2 distribution')
mu2 = createdistribution(n2,5*pi/6,0.2)
n3 = 80;
disp('Mu3 distribution')
mu3 = createdistribution(n3,3*pi/2,0.2)
allpoints = vertcat(mu1, mu2, mu3)
sampleSizes = [n1 n2 n3]';
baryCenterPoints=barycenterndistributions(allpoints, sampleSizes)

sampleShapes = ['s';'d';'o'];
sampleColors = ['y';'g';'b'];
h2 = subplot(1,3,2);
titletext= "Three distributions";
plotresults(sampleSizes, sampleShapes, sampleColors,...
    allpoints, baryCenterPoints, titletext)

n1 = 60;
disp('Mu1 distribution')
mu1 = createdistribution(n1,0,0.2)
n2 = 80;
disp('Mu2 distribution')
mu2 = createdistribution(n2,pi/4,0.2)
n3 = 80;
disp('Mu3 distribution')
mu3 = createdistribution(n3,3*pi/4,0.2)
n4 = 80;
disp('Mu4 distribution')
mu4 = createdistribution(n4,5*pi/4,0.2)
allpoints = vertcat(mu1, mu2, mu3, mu4)

sampleSizes = [n1 n2 n3 n4]';
baryCenterPoints=barycenterndistributions(allpoints, sampleSizes)
sampleShapes = ['s';'d';'o';'+';'*'];
sampleColors = ['y';'m';'g';'b'];
h2 = subplot(1,3,3);
titletext= "Four distributions";
plotresults(sampleSizes, sampleShapes, sampleColors,...
    allpoints, baryCenterPoints, titletext)

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'barycenter','-dpdf','-r0')


function xy=createdistribution(n, theta0, sigma)
theta = theta0 + sigma.*randn(n,1);
xy = zeros(n,2);
xy(1:n,1) = cos(theta);
xy(1:n,2) = sin(theta);
end

function plotresults(sampleSizes, sampleShapes, sampleColors,...
    allpoints, baryCenter, titletext)
nssz = size(sampleSizes,1);
n= 0;
nprev = 0;
for i=1:nssz
    nsampleSizes = sampleSizes(i);
    n = n + nsampleSizes;
    pl = plotpoints(allpoints(1+nprev:n,1),allpoints(1+nprev:n,2),sampleShapes(i), 7, sampleColors(i));
    set(get(get(pl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    hold on
    nprev = n;
end
plotpoints(baryCenter(:,1), baryCenter(:,2), 'p', 10, 'r');
axis square
% lgd=legend('Barycenter','Location','southeast');
% lgd.FontSize = 10;
% lgd.TextColor = 'black';
xlabel('X')
ylabel('Y')
title(titletext, 'FontSize', 14);
hold off;
end


function pl=plotpoints(x, y, marker, size, col)
myplot = @(x,y,col)plot(x,y,marker,...
    'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', col,...
    'MarkerSize', size,...
    'LineWidth', 2);
pl = myplot(x, y, col);
end

