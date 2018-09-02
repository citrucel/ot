function [minT, maxT] =computeminmaxhourly(sampleSizes, temperatures)

sz = length(sampleSizes);

minT = zeros(sz,4);
maxT = zeros(sz,4);
j= 1;
for i=4:4:16
    [min, max] = computeminmaxhourlywithhour(sampleSizes, temperatures, i);
    minT(:,j) = min;
    maxT(:,j) = max;
    j = j + 1;
end

end

function [minT, maxT] = computeminmaxhourlywithhour(sampleSizes, temperatures, hour)
[m,~] = size(sampleSizes);
minT = zeros(m,1);
maxT = zeros(m,1);
sz = 0;
for i=1:m
    nsz = sz + sampleSizes(i);
    temperatureMonths = temperatures(sz+1:nsz, hour);
    temperatureMonths = temperatureMonths(:);
    minT(i)= min(temperatureMonths);
    maxT(i)= max(temperatureMonths);
    sz = nsz;
end
end