function stats =find_statistical_measures(var)

% find mean, median, 25th, 75th percentile and standard deviation for both 
% input variables;
nvar = size(var,1);
for i = 1:nvar
    stats(i).mean = mean(var(i,:),'omitnan');
    stats(i).median = median(var(i,:), 'omitnan');
    stats(i).InterQuartile = prctile(var(i,:), [25, 75]);
    IQR = iqr(var(i,:));
    stats(i).upperWhisker = stats(i).InterQuartile(2)+1.5*IQR;
    stats(i).lowerWhisker = stats(i).InterQuartile(1)-1.5*IQR;
    stats(i).stdv = std(var(i,:), 1,'omitnan');
end


return