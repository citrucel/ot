clc;
clear;

n0 = 10
n1 = 10
[xy,  xypairs, expxymu] = optimaltransport(n0, n1, 0);

h = figure;
h1 = subplot(1,4,1);
plotresults(n0, xy, xypairs, expxymu);
h2 = subplot(1,4,2);
n0 = 60
n1 = 80
[xy,  xypairs, expxymu] = optimaltransport(n0, n1, 0);
plotresults(n0, xy, xypairs, expxymu);
h3 = subplot(1,4,3);
n0 = 120
n1 = 140
[xy,  xypairs, expxymu] = optimaltransport(n0, n1, 0);
plotresults(n0, xy, xypairs, expxymu);

h4 = subplot(1,4,4);
[xy,  xypairs, expxymu] = optimaltransport(n0, n1, 1);
plotresults(n0, xy, xypairs, expxymu);

% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,'twodist','-dpdf','-r0')

function [xy,  xypairs, expxymu] = optimaltransport(n0, n1, filterParameter)
n = n0 + n1;
xy = zeros(n,2);
rho = randn(1,n0)*1;
xy(1:n0,1) = cos(rho);
xy(1:n0,2) = sin(rho);
disp('Rho distribution')
xy(1:n0,:)
mu = pi/2 + randn(1,n1)*0.5;
xy(n0+1:n,1) = cos(mu);
xy(n0+1:n,2) = sin(mu);
disp('Mu distribution')
xy(n0+1:n,:)


dmat = computedistances(xy);
disp('Initial distance matrix')
dmat

%plotdistmatrix(xy,dmat)
D=dmat(1:n0, n0+1:n);
%D=dmat;
sort(D(:));
xypairs=findclosestpairofpoints(D);

mindmat=mindist(dmat, filterParameter);
disp('Minimal distance matrix')
mindmat

sqMinDistMat= (mindmat.^2)./2;
disp('Squared minimal distance matrix')
sqMinDistMat

disp('Minimal distance cost matrix')
minCost = sqMinDistMat(1:n0, n0+1:n);
minCost
if isOneInfinite(minCost)
    return
end

p0 = ones(n0,1) * 1/n0;
p1 = ones(n1,1) * 1/n1;
pimatrix = computewithlinprog(minCost, p0, p1);
expxymu = (n0 .* pimatrix) * xy(n0+1:n, :);
end

function [minxy]=findclosestpairofpoints(A)
n = size(A,1);
m = size(A,2);
nn = min(n,m);
mm = max(n,m);
B=reshape(A,nn,mm);
minxy = zeros(nn,2);
C=B(:,:);
for i=1:nn
    [minMatrix, ~] = min(min(C));
    [row,col] = find(B==minMatrix);
    minxy(i,1)=row;
    minxy(i,2)=col;
    [row,col] = find(C==minMatrix);
    C(:,col)=[];
    C(row,:)=[];
end
end

function plotresults(n0, xy, xypairs, xymu)
n = size(xy,1);
pl1 = plotpoints(xy(1:n0,1),xy(1:n0,2),'r');
set(get(get(pl1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
hold on
pl2 = plotpoints(xy(n0+1:n,1),xy(n0+1:n,2),'g');
set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
hold on
pl3 = plotpoints(xymu(:,1), xymu(:,2),'b');
set(get(get(pl3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
for k=1:n0
    xi = xypairs(k,1);
    yi = xypairs(k,2);
    plot([xy(xi,1) xy(n0+yi,1)], [xy(xi,2) xy(n0+yi,2)], 'b', 'LineStyle', '--', 'LineWidth', 0.5);
    plot([xy(k,1) xymu(k,1)], [xy(k,2) xymu(k,2)], 'g' , 'LineStyle', '--', 'LineWidth', 0.5);
end
axis square
lgd=legend('min d(X,Y)', 'T(X)','Location','southwest');
lgd.FontSize = 10;
lgd.TextColor = 'black';
txt = texlabel(['|X|=' num2str(n0) ' and |Y|=' num2str(n-n0)]); 
title(txt, 'FontSize', 14)
end


function pl = plotpoints(x, y, col)
myplot = @(x,y,col)plot(x,y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);
pl = myplot(x, y, col);
end

function plotdistmatrix(xy, dmat)
figure; plot(xy(:,1),xy(:,2),'.');
n = size(xy,1);
for i=1:n, text(xy(i,1),xy(i,2),[' ' num2str(i)]); end
figure; imagesc(dmat); colorbar
end





