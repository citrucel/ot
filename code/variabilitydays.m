function maxMinT = variabilitydays(month)
maxMinT = -Inf;
[ndays,~] = size(month);
for i=1:ndays
    minT = min(month(i,:));
    maxT = max(month(i,:));
    d = maxT - minT;
    if d >= maxMinT
        maxMinT = d;
    end
end
end