function hourlyData = computebyhour(temperatures)
[~,n] = size(temperatures);
hourlyData = zeros(3,n);
for i=1:n
    hourlyTemperatures = temperatures(:,i);
    hourlyData(1,i)=min(hourlyTemperatures);
    hourlyData(2,i)=median(hourlyTemperatures);
    hourlyData(3,i)=max(hourlyTemperatures);
end
end