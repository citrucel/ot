function [coldestDay, warmestDay,...
    maxMinMonth, meanMonth, stdMonth,...
    hourlyData] = computebarycenterstats(baryCenterPoints)
[coldestDay,warmestDay] = coldestwarmestdays(baryCenterPoints);
maxMinMonth = variabilitydays(baryCenterPoints);
[meanMonth, stdMonth] = computemeanstd(baryCenterPoints);
hourlyData = computebyhour(baryCenterPoints);
end

function [meanMonth, stdMonth] = computemeanstd(month)
temperatureMonths = month(:);
meanMonth = mean(temperatureMonths);
stdMonth =std(temperatureMonths);
end
