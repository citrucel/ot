

clc;
clear;
load('AdjustedPreparedCRNH02032016GANewton8W2.mat', 'data')
seasons = data.season;
idx = double(seasons) ;
pointMarkers = arrayfun(@Season.getSeasonMarker,idx);
% y = randi(10, 365, 1);
% x = randi(10, 365, 1);
% figure
% scatter (x(:,1), y(:,1), 'MarkerFaceColor', [102, 215, 209]./255)
% hold on

x = linspace(0,3*pi,365)';
y = cos(x) + rand(1,365)';
sz = 25;
c = linspace(1,10,length(x));
[winter,spring,summer,fall]=partitionBySeason(x, y, idx);
 scatter(winter(:,1),winter(:,2),sz,'Marker','o')
 
% for i=1:size(winter,1)
%     scatter(winter(i,1),winter(i,2),sz,'Marker','o')
% end

function [winter,spring,summer,fall]=partitionBySeason(x, y, seasons)
n = size(x,1);
winter = zeros(n,2);
spring = zeros(n,2);
summer = zeros(n,2);
fall = zeros(n,2);
j = 1;
for i=1:n
    % Winter
    if seasons(i) == 1
        winter(j,1) = x(i);
        winter(j,2) = y(i);
        % Spring
    elseif seasons(i) == 2
        spring(j,1) = x(i);
        spring(j,2) = y(i);
        % Summer
    elseif seasons(i) == 3
        summer(j,1) = x(i);
        summer(j,2) = y(i);
        % Fall
    elseif seasons(i) == 4
        fall(j,1) = x(i);
        fall(j,2) = y(i);
    end
    j = j +1;
end


end
% xy(:,1)  = randi(10,365,1);
% xy(:,2)  = randi(10,365,1);
% xy(:,1) = linspace(0,3*pi,200);
% xy(:,2) = linspace(0,3*pi,200);
% xy = zeros(365,2);
% rho = randn(1,365)*1;
% xy(:,1) = cos(rho);
% xy(:,2) = sin(rho);
% c=[102, 215, 209;253,174,97;215,25,28]/255;
% m=['o','d','s'];
% figure, hold on
% for i=1:size(x,2)
%   scatter (x(:,i), y(:,i), 'MarkerFaceColor', c(i,:),'Marker',pointMarkers(i))
% end

% figure
% for i=1:365
%     %pointMarkers(i)
%     plot(xy(i,1), xy(i,2), 'r');
%     %plotpoints(xy(i,1), xy(i,2), pointMarkers(i), 0.5, 'r');
%     %plot(xy(i,1), xy(i,2))
%     %scatter(xy(i,1),xy(i,2), 'Marker',pointMarkers(i))
% end
% hold on

% sampleSizes = data.sampleSizes;
% global allpoints;
% allpoints = data.allpoints;
% seasons = data.season;
%
% global sqMinDistMat1;
% sqMinDistMat1 = computesquareddistances(allpoints, true, 0);
%
% %Y = tsne(allpoints,'Algorithm','barneshut','Distance',@distancethrehsold1);
% Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold1);
% subplot(3,1,1)
% gscatter(Y(:,1),Y(:,2),seasons')
%
% global sqMinDistMat2;
% sqMinDistMat2 = computesquareddistances(allpoints, true, 1);
% Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold2);
% subplot(3,1,2)
% gscatter(Y(:,1),Y(:,2),seasons')
%
% global sqMinDistMat3;
% sqMinDistMat3 = computesquareddistances(allpoints, true, 2);
% Y = tsne(allpoints,'Algorithm','exact','Distance',@distancethrehsold3);
% subplot(3,1,3)
% gscatter(Y(:,1),Y(:,2),seasons')

%
% load fisheriris
%
% rng('default') % for reproducibility
% Y = tsne(meas,'Algorithm','exact','Distance','mahalanobis');
% subplot(2,2,1)
% gscatter(Y(:,1),Y(:,2),species)
% title('Mahalanobis')
%
% rng('default') % for fair comparison
% Y = tsne(meas,'Algorithm','exact','Distance','cosine');
% subplot(2,2,2)
% gscatter(Y(:,1),Y(:,2),species)
% title('Cosine')
%
% rng('default') % for fair comparison
% Y = tsne(meas,'Algorithm','exact','Distance','chebychev');
% subplot(2,2,3)
% gscatter(Y(:,1),Y(:,2),species)
% title('Chebychev')
%
% rng('default') % for fair comparison
% Y = tsne(meas,'Algorithm','exact','Distance','euclidean');
% subplot(2,2,4)
% gscatter(Y(:,1),Y(:,2),species)
% title('Euclidean')

% R = 10; x_c = 5; y_c = 8;
% thetas = 0:pi/64:pi;
% xs = x_c + R*cos(thetas);
% ys = y_c + R*sin(thetas);
%
% % Now add some random noise to make the problem a bit more challenging:
% mult = 0.5;
% xs = xs+mult*randn(size(xs));
% ys = ys+mult*randn(size(ys));
%
% xy
% figure
% plot(xs,ys,'b.')

% n0=1000
% th = linspace( 0, 2*pi, n0 ); %// N samples
% noise = rand( 1, n0 ) * .1; %// random noise in range [0..0.1]
% r = 1+noise; %// add noise to r=1
% figure;
% plot( r.*cos(th), r.*sin(th) ); title('noise on circle');

% n0=2000
% xy = createdistribution(n0, 0, 10);
% plot(xy(:,1),xy(:,2),'b.')
% axis equal


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

function xyn=createdistribution(n, theta0, sigma)
theta = theta0 + sigma.*randn(1,n);
xy = zeros(n,2);
xy(1:n,1) = cos(theta);
xy(1:n,2) = sin(theta);
noise = rand( 1, n) * .1;
r = 1 + noise;
xyn = zeros(n,2);
xyn(1:n,1) = r .* cos(theta);
xyn(1:n,2) = r .* sin(theta);
end



% xy = createdistribution(2000, 0, 1);
% %plot(xy(:,1), xy(:,2),'p', 10, 'r');
% plotpoints(xy(:,1), xy(:,2), 'p', 10, 'r');

% load('PreparedCRNH02032016GANewton8W.mat', 'data')
% sampleSizes = data.sampleSizes;
% allpoints = data.allpoints;
% allpoints=allpoints(1:5,:);

%allpoints = [0 2 3; 2 0 4; 3 4 0]
% distanceMat = computedistances(allpoints)
% distanceMat2 = computedistances2(allpoints)
%A=filtermatrix(allpoints)

% function xy=createdistribution(n, theta0, sigma)
% theta = theta0 + sigma.*randn(n,1);
% xy = zeros(n,2);
% xy(1:n,1) = cos(theta);
% xy(1:n,2) = sin(theta);
% end

function pl=plotpoints(x, y, marker, size, col)
myplot = @(x,y,col)plot(x,y,marker,...
    'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', col,...
    'MarkerSize', size,...
    'LineWidth', 2);
pl = myplot(x, y, col);
end

function new=filtermatrix(A)
[n,m]=size(A);
new = inf(n, m);
idx = 1==eye(size(A));
new(idx) = 0;
B = A;
B(idx) = inf;
[minv,idx] = min(B, [], 1);
for i=1:n
    j = idx(i);
    new(i, j) = minv(i);
    new(j,i) = new(i,j);
end
end

function distanceMat = computedistances(allpoints)
n = size(allpoints,1);
distanceMat = zeros(n,n);
for i=1:n
    for j=i+1:n
        distanceMat(i,j) = pdist2(allpoints(i,:), allpoints(j,:),'euclidean');
        distanceMat(j,i) = distanceMat(i,j);
    end
end
end

function distanceMat = computedistances2(allpoints)
[n,m] = size(allpoints);
distanceMat = zeros(n,n);
for i=1:n
    for j=i+1:n
        %             d  = sqrt(sum((allpoints(i,k) - allpoints(j,k)) .^ 2));
        d  = norm(allpoints(i,:) - allpoints(j,:));
        distanceMat(i,j) = d;
        distanceMat(j,i) = d;
    end
end
end

% figure
% h = subplot(4,4,[1 2 5 6])
% bar(normrnd(0,1,100,1))
% subplot(4,4,[3,4,7,8])
% histogram(normrnd(0,1,100,1))
% subplot(4,4,[9,10,13,14])
% bar(normrnd(0,1,100,1))
% subplot(4,4,[11,12,15,16])
% histogram(normrnd(0,1,100,1))

% load('NewtonT.mat', 'saveStruct')
% sampleSizes = saveStruct.sampleSizes;
% allpoints = saveStruct.allpoints;
% barycenter = saveStruct.baryCenter;
% [m,n] = size(allpoints)
% [m2,n2] = size(barycenter)
% dayst = zeros(m + m2,n)
% dayst(1:m,1:n) = allpoints(:,:)
% dayst(m+1:m+m2,:)=barycenter(1:m2,:)
% xlswrite('dayst.xls',dayst)
% M = csvread('cleaneddayst.csv')

% doMap=false;
% PlacesYear;

% XX: matrix with data from year read.

% nst=length(places);

% year_reg=[31 28 31 30 31 30 31 31 30 31 30 31];
% year_leap=0*year_reg; year_leap(2)=1;
% c_reg=[0 cumsum(year_reg)];
% c_leap=[0 cumsum(year_leap)];
% leap_years=[2000 2004 2008 2012 2016];
% Years=[2000:2017];
% ndays=365+ismember(Years,leap_years);
% c_days=[0 cumsum(ndays)];

%for jst=45:nst,

% file_name='CRNH0203-2016-GA_Newton_8_W.txt';
%
%     cd('../../CRNH0203-201508171551')
%     cd 2016
%     file_name=strcat('CRNH0203-2016-',places{jst}{5},'.txt');
% clearvars VarName*
% uiimport(file_name);

% import(file_name);

% T = readtable('CRNH0203-2016-GA_Newton_8_W.txt','Delimiter',' ');


% pause

%     cd('../../Current/Boulder')
%
% m=length(VarName2);
%
% XX=zeros(m,5);
%
% D=VarName2; % Date.
% T=VarName3; % Time
%
% year=floor(D/10000);
% month=floor((D-year*10000)/100);
% day=floor((D-year*10000-month*100));
%
% ItsLeap=ismember(year,leap_years);
%
% t=(c_days(year-1999)'+c_reg(month)'+ItsLeap.*c_leap(month)'+day)*24+T/100; % Time in hours from Jan 1st 2000 at 0am.
%
% XX(:,1)=t/24;      % Time in days from Jan. 1st 2000 at 0am.
%
% XX(:,2)=VarName10; % Temperature;
%
% XX(:,3)=VarName13; % Precipitation;
%
% XX(:,4)=VarName14; % Solar radiation;
%
% XX(:,5)=VarName15; % Flag for it (0: good, 3: bad);
%
% file_name=strcat('Data.mat');
% %     file_name=strcat('Data',places{jst}{1},'.mat');
%
% load(file_name);
%
% t_old=max(XXX(:,1));
%
% XXX=[XXX; XX(XX(:,1)>t_old,:)];
%
% save(file_name,'XXX')

%end



%
% v = zeros(3,1)
% for i=1:3
%  v(i)= rand;
% end
% v
%
% A = rand(4,5)
% B = rand(4,7)
% C=cell(1);
% C{1,1} = A;
% C{1,2} = B;
% C
% AT = horzcat(C{:})

% figure(1)
% surf(peaks);
% colormap(winter);
% title('FIGURE 1A', 'FontSize', 12, 'fontweight', 'bold')
%
% figure(2)
% surf(peaks);
% colormap(autumn);
% title('FIGURE 1B', 'FontSize', 12, 'fontweight', 'bold')
%
% figure(3)
% surf(peaks);
% colormap(spring);
% title('FIGURE 1C', 'FontSize', 12, 'fontweight', 'bold')
%
% export_fig FIGURE_1.tiff -m3 -q101 -nocrop



% h1 = subplot(1,3,1);
% surf(peaks);
% colormap(h1, winter);
% axis square;
% title('FIGURE 1A', 'FontSize', 12, 'fontweight', 'bold')
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%
% h2 = subplot(1,3,2);
% surf(peaks);
% axis square;
% colormap(h2, autumn);
% title('FIGURE 1B', 'FontSize', 12, 'fontweight', 'bold')
%
% h3 = subplot(1,3,3);
% surf(peaks);
% colormap(spring);
% axis square;
% title('FIGURE 1C', 'FontSize', 12, 'fontweight', 'bold')




% f = [66.44 13.8 8.4];
% A = [0.4 4.3 1.85; 0.377 0.43 0.5];
% b = [32 5];
% Aeq = [1 1 1];
% beq = [9];
% lb = [1 1 1];
% ub = [9 9 9];
%
% [XYZ,fval]=linprog(f,A,b,Aeq,beq,lb,ub)
% f = rand(6, 1)
% A=[1 1 0 0 0 0; ...
%     0 0 0 1 1 1;...
%     1 0 0 1 0 0;
%     0 1 0 0 1 0;
%     0 0 1 0 0 1];
% A=[    1     0     1     0     1     0;...
%      0     1     0     1     0     1;...
%      1     1     0     0     0     0;...
%      0     0     1     1     0     0;...
%      0     0     0     0     1     1];
%
% %     1 1 1 1 1 1];
% b = [1/2 1/2 1/3 1/3 1/3];
% Aeq = [1 1 1 1 1 1];
% beq = [1];
% % lb = [-Inf -Inf -Inf -Inf -Inf -Inf];
% % ub = [Inf Inf Inf Inf Inf Inf];
% lb = [0 0 0 0 0 0];
% ub = [1 1 1 1 1 1];
% [XYZ,fval]=linprog(f,[],[],A,b,lb,ub)
% reshape(XYZ,2,3)
%,lb,ub)

%
% A = [2 6 5 3; 7 -2 4 6; 1 0 -9 9];
% [v,d] = min(A,[],2);
% sort(A(:))
% minxy=findclosestpairofpoints(A)
%
%
% function [minxy]=findclosestpairofpoints(A)
%     n = size(A,1);
%     m = size(A,2);
%     nn = min(n,m);
%     minxy = zeros(nn,2);
%     B=A(:,:);
%     for i=1:nn
%         [minMatrix, ~] = min(min(B));
%         [row,col] = find(A==minMatrix);
%         minxy(i,1)=row;
%         minxy(i,2)=col;
%         [row,col] = find(B==minMatrix);
%         B(:,col)=[];
%         B(row,:)=[];
%     end
% end

% n0 = 3
% n1 = 4
% n = n0 + n1;
% xy = zeros(n,2);
% rho = randn(1,n0)*1;
% xy(1:n0,1) = cos(rho);
% xy(1:n0,2) = sin(rho);
% disp('Rho distribution')
% xy(1:n0,:)
% mu = pi/2 + randn(1,n1)*0.5;
% xy(n0+1:n,1) = cos(mu);
% xy(n0+1:n,2) = sin(mu);
% disp('Mu distribution')
% xy(n0+1:n,:)
%
% dmat = computedist(xy(:,1),xy(:,2));
% disp('Initial distance matrix')
% dmat
%
% mindmat=mindist(dmat)
%(mindmat.^2)./2

% disp('Minimal distance matrix')
% niter = 20;
% A= zeros(4,4);
% A(:,:) = inf;
% A(1,1) = 0;
% A(1,3) = -2;
% A(2,1) = 4;
% A(2,2) = 0;
% A(2,3) = 3;
% A(3,3) = 0;
% A(3,4) = 2;
% A(4,2) = -1;
% A(4,4) = 0;
%A
% B = mindist(A, 5)

% A= rand(4,4);
% A(1,1)=0;
% A(2,2)=0;
% A(3,3)=0;
% A(4,4)=0;
% for i=1:4
%     for j=i+1:4
%         A(j,i)=A(i,j);
%     end
% end

% A= zeros(7,7);
%
%  A=    [         0    0.0746    0.0000    1.7411    0.6862    1.9946    0.7760;
%     0.0746         0    0.0734    1.4314    0.3498    1.8810    0.4234;
%     0.0000    0.0734         0    1.7390    0.6832    1.9942    0.7729;
%     1.7411    1.4314    1.7390         0    0.5950    0.1931    0.5117;
%     0.6862    0.3498    0.6832    0.5950         0    1.2133    0.0043;
%     1.9946    1.8810    1.9942    0.1931    1.2133         0    1.1214;
%     0.7760    0.4234    0.7729    0.5117    0.0043    1.1214         0]
% B = mindist(A)



% A= zeros(4,4);
% A(:,:) = inf;
% A(1,3) = -2;
% A(3,4) = 2;
% A(4,2) = -1;
% A
% B = computemindist(A)

% n0 = 2000
% n1 = 2000
% n = n0 + n1;
% xy = zeros(n,2);
% rho = randn(1,n0)*1;
% xy(1:n0,1) = cos(rho);
% xy(1:n0,2) = sin(rho);
% disp('Rho distribution')
% xy(1:n0,:)
% mu = pi/2 + randn(1,n1)*0.5;
% xy(n0+1:n,1) = cos(mu);
% xy(n0+1:n,2) = sin(mu);
% disp('Mu distribution')
% xy(n0+1:n,:)


%
% figure
% plotpoints(xy(1:n0,1),xy(1:n0,2),'r')
% hold on
% plotpoints(xy(n0+1:n,1),xy(n0+1:n,2),'g')
% hold on

function A = computedist(x, y)
np = size(x,1);
A = zeros(np,np);
for i=1:np
    for j=i+1:np
        A(i,j)=((x(i)-x(j))^2+(y(i)-y(j))^2)^(1/2);
        A(j,i) = A(i,j);
    end
end
end

function B = computemindist(A)
n = size(A, 1);
B = A(:,:);
for k=1:n
    for i=1:n
        for j=1:n
            s = B(i,k) + B(k,j);
            if B(i,j) > s
                B(i,j) = s;
            end
        end
    end
end
end

