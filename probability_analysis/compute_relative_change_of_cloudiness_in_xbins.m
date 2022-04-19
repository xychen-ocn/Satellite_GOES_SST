function [RCC, RCC_domain_mean, xedges] = compute_relative_change_of_cloudiness_in_xbins(x, xedges_prc, daily_CF, mean_CF, LON, LAT)
% This script will include a function to compute the relative change in
% temporal cloud fraction (or cloudiness) for each bin of effective SST
% gradient or SST laplacian. 
% follow the same idea as in Des.. et al. (2021).

% input: x (SST gradient or laplacian)
%        xedges_prc (percentile)
%        cldmasks (?)
%        mean_cldfreq;
% update this function on Apr 19, 2022:

% check:
if size(x)==size(daily_CF)

nedges = length(xedges_prc);

xvec = reshape(x, 1,[]);
xedges = prctile(xvec, xedges_prc);

%CFvec = reshape(daily_CF, 1,[]);
%LONs = repmat(LON, 1,1, size(x,3));
%LATs = repmat(LAT, 1,1, size(x,3));
if length(size(mean_CF))==2
    % mean is only temporal:
    spatioTemp_mean = False;
    mean_CF_3D = repmat(mean_CF, 1,1,size(x,3));
elseif length(size(mean_CF))==1
    % mean is both spatial and temporal:
    spatioTemp_mean = True;
    
end
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
        if ~spatioTemp_mean
            CF_change = (daily_CF(inds) - mean_CF_3D(inds))./mean_CF_3D(inds) * 100;
            CF_change_domain_mean = (daily_CF(inds) - mean(mean_CF(:)))./mean(mean_CF(:)) * 100;
            
        else
            CF_change = (daily_CF(inds) - mean_CF)./mean_CF * 100;
            CF_change_domain_mean = CF_change;
        end
        
        RCC(id) = mean(CF_change,'omitnan');
        RCC_domain_mean(id) = mean(CF_change_domain_mean,'omitnan');

           
    else
        disp('no data found for this bin');
        RCC(id) = NaN;
        RCC_domain_mean = NaN;
    end
    
end

else
    disp('inconsistent size between x and cloudiness variable');
    return
    
end
    




end
