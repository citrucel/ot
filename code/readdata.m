clc
clear
load('AdjustedCRNH02032016GANewton8W.mat')
DataTable=CRNH02032016GANewton8W;
clear AdjustedCRNH02032016GANewton8W

ts = length(DataTable.LocalDay);
sampleSizes = zeros(12,1);
[prevMonth,prevDay] = convertDate(DataTable.LocalDay(1));
allpoints = zeros(365,24);
seasonlabels{365} = [];
% seasonlabels = zeros(1,365);
temperatures = zeros(1,24);
j=1;
k=1;
sn = 1;
nd = 1;
for i=1:ts
    D = DataTable.LocalDay(i);
    %     if D > 20170000
    %         T = DataTable.T(i);
    %         temperatures(j)=T;
    %         allpoints(nd,1:24) = temperatures(:);
    %         season = Season.getSeason(currMonth);
    %         seasonlabels{nd} = string(season);
    %         sampleSizes(k)=sn;
    %         break;
    %     end
    [currMonth,currDay] = convertDate(D);
    if currDay == prevDay
        T = DataTable.T(i);
        temperatures(j)=T;
        j = j + 1;
    else
        season = Season.getSeason(currMonth);
        seasonlabels{nd} = season;
        allpoints(nd,1:24) = temperatures(:);
        nd = nd + 1;
        temperatures = zeros(1,24);
        temperatures(1) = DataTable.T(i);
        j = 2;
        prevDay = currDay;
        if currMonth == prevMonth
            sn = sn + 1;
        else
            prevMonth = currMonth;
            sampleSizes(k)=sn;
            k = k + 1;
            sn = 1;
        end
    end
end
%allpoints(allpoints==-9999)=0;
allpoints(nd,1:24) = temperatures(:);
% season = Season.getSeason(currMonth);
% seasonlabels(nd) = season;
season = Season.getSeason(currMonth);
seasonlabels{nd} = season;
sampleSizes(k)=sn;
data.sampleSizes = sampleSizes;
data.allpoints = allpoints;
seasonCategories=categorical(seasonlabels,{'Winter','Spring','Summer','Fall  '},'Ordinal',true);
data.season = seasonCategories;
save('AdjustedPreparedCRNH02032016GANewton8W2.mat','data')

function [month, day]=convertDate(D)
md =  D - 20160000;
month = fix(md/100);
day = rem(md,100);
end



