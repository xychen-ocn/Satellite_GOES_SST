function [RCC, RCC_domain_mean, xedges] = compute_relative_change_of_cloudiness_in_xbins(x, xedges_prc, daily_CF, mean_cldfreq, LON, LAT)
% This script will include a function to compute the relative change in
% temporal cloud fraction (or cloudiness) for each bin of effective SST
% gradient or SST laplacian. 
% follow the same idea as in Des.. et al. (2021).

% input: x (SST gradient or laplacian)
%        xedges_prc (percentile)
%        cldmasks (?)
%        mean_cldfreq;

nedges = length(xedges_prc);

xvec = reshape(x, 1,[]);
xedges = prctile(xvec, xedges_prc);

%CFvec = reshape(daily_CF, 1,[]);
%LONs = repmat(LON, 1,1, size(x,3));
%LATs = repmat(LAT, 1,1, size(x,3));
mean_cldfreq_3D = repmat(mean_cldfreq, 1,1,size(x,3));
%% 1. discretize x data into different bins 
GrpIDs = discretize(x,xedges);
RCC = zeros(1, nedges-1);
RCC_domain_mean = RCC;
for id = 1:nedges - 1
    inds = find(GrpIDs==id);   % locations of x;
    if ~isempty(inds)
        % use these indices to compute cloudiness change:
        %xlocs = LONs(inds);
        %ylocs = LATs(inds);
        CF_change = (daily_CF(inds) - mean_cldfreq_3D(inds))./mean_cldfreq_3D(inds) * 100;
        CF_change_domain_mean = (daily_CF(inds) - mean(mean_cldfreq(:)))./mean(mean_cldfreq(:)) * 100;

        RCC(id) = mean(CF_change,'omitnan');
        RCC_domain_mean(id) = mean(CF_change_domain_mean,'omitnan');
        
    else
        disp('no data found for this bin');
        RCC(id) = NaN;
        RCC_domain_mean = NaN;
    end
    
end
    




end
