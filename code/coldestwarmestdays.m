function [coldestDay,warmestDay] = coldestwarmestdays(month)
coldestDay = +Inf;
warmestDay = -Inf;
[ndays,~] = size(month);
for i=1:ndays
    meanDayT = mean(month(i,:));
    if meanDayT <= coldestDay
        coldestDay = meanDayT;
    end
    if meanDayT >= warmestDay
        warmestDay = meanDayT;
    end
end
end