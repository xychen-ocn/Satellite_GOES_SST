function null = bootstrap_nullhypothesis_level_for_relative_cloudiness_change( x, xprc, daily_CF, mean_cldfreq, varargin)

switch nargin
    case 4
        nrep = 1000;
    case 5
        nrep = varargin{1};
end

if xprc > 1
    ndata = numel(x)*xprc*0.01;
else
    ndata = numel(x)*xprc;
end

% sample 1000 times and take out ndata from x randomaly.
mean_cldfreq_3D = repmat(mean_cldfreq, 1,1,size(x,3));

for it = 1:nrep
    
    ridx = randperm(numel(x), ndata);
    
    %x_samples = x(ridx);
    daily_CF_samples = daily_CF(ridx);
    
    CF_change = (daily_CF_samples - mean_cldfreq_3D(ridx))./mean_cldfreq_3D(ridx)*100;
    
     
    mean_CF_change(it) = mean(CF_change(:), 'omitnan');
    stdv_CF_change(it) = std(CF_change(:), 1, 'omitnan');
    
    CF_change_samples(:,it) = CF_change;
    
    if mod(round(it/nrep*100), 10)==0
        disp(['finished ' num2str(round(it/nrep*100)) '% sampling...'])
    end
    
end

% figure(10);
% subplot(3,1,1);
% histogram(CF_change_samples(:));
% 
% subplot(3,1,2);
% histogram(mean_CF_change(:));
% 
% subplot(3,1,3);
% histogram(stdv_CF_change(:));

null.mean = mean(mean_CF_change, 'omitnan');
null.stdv = std(mean_CF_change, 1, 'omitnan');
null.p95 = prctile(mean_CF_change, 95);
null.p05 = prctile(mean_CF_change, 5);

    


return