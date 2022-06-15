function cloudInfo = compute_cloudfrac_over_features_from_L3C_cloudmasks(cloudfrac_data,t, search_box, blobsIn, cutoff_thres,varargin)
% purpose: this script is used to compute cloudiness from L3C cloudmasks,
% in the region of interest. as well as computing the cloudiness
% enhancement from a reference level. 
% Inputs: 
%     t: time of request;
%     search_box: area to compute cloudiness from L3C cloudmask; 
%  - need to deal with the fact that L3C data is about 2.5x finer than the
%  L4 data; (how do we use the cloud mask information for a coarser grid?
%  In this function, the cloudiness is computed at the native resolution (2km).)
%
%
  %%%%%%%%%%%% Priority %%%%%%%%%%%%%%%
  % I wanted to have the 2-month cloud frequency at the location of each
  % eddies. need to make this (2-month cloud frequency data stored for access.)
  % will make the 2-month mean cloud frequency data stored for access.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin 
    case 4
        masktype = 'circle';
    case 5
        masktype = 'square';
end
  
% load in cloud fraction data:
 %cloudfrac_data = './GOES16_L3C_cloudfrac/cloudfrac_5kmres_calculated_from_2km_cloudmask_8-18N_62W-45W.mat';
 
 load(cloudfrac_data); % contains LON_cldfrac, LAT_cldfrac, cloudiness_ds;
 
%  if strcmpi(refstate, 'JanFebMean')
%     tmp = cat(3,cloudiness_ds.cldfrac_hrly);
%     ref_CF = mean(tmp, 3, 'omitnan');
%     
% elseif strcmpi(refstate, '3d_movmean')
%     %ref_CF = CF_refdata.CF_movmean_2monthsAve.wdsz_3d;
%     % need to do the moving mean manually.
% elseif strcmpi(refstate, '5d_movmean')
%     %ref_CF = CF_refdata.CF_movmean_2monthsAve.wdsz_5d;
%     % need to do the moving mean manually.
% else
%     disp('invalid refstate, please enter either *JanFebMean*, *3d_movmean*, or *5d_movmean*. For more options, please generate the reference state first!')
%     return
% end

% load in L3C daily dataset;
tidx = find(time_cldfrac==t);  % find out the time index to take the cloudfrac results out.

% daily cloud fraction:
nhrly_maps = size(cloudiness_ds(tidx).cldfrac_hrly,3);
if nhrly_maps>=12
    cldfrac_today = cloudiness_ds(tidx).cldfrac_dailymean;    % this is already the best_quality
else
    cldfrac_today = nan(size(LON_cldfrac));                   % exclude Feb 26 for example.
end

% highpass filter cloud fraction with different filtering days:
NT = length(time_cldfrac);
cldfrac_2mon_d = reshape([cloudiness_ds.cldfrac_dailymean], size(LON_cldfrac,1), [], NT);

% cutoff day will be inputed.
% -- the following is temporal filtering (tested, it does not work as
% expected.)
%daily_cldfrac_highpassed = highpass_filter_large_scale_cloudiness(cldfrac_2mon_d, cutoff_day); 
for it = 1:length(time_cldfrac)
    daily_cldfrac_LS(:,:,it) = estimate_LSG_SST(LON_cldfrac(1,:), LAT_cldfrac(:,1), cldfrac_2mon_d(:,:,it), 'method','spectrum','CutoffScale', cutoff_thres, 'checkflag',false);
end
cldfrac_highpassed_today = cldfrac_2mon_d(:,:,tidx)- daily_cldfrac_LS(:,:,tidx);
cldfracbg_today = daily_cldfrac_LS(:,:,tidx);    % daily background cldfrac field.
refCF = mean(daily_cldfrac_LS,3);       % 2 month background cldfrac.
%% QC data:
% throw away maps that is 100% garbage; already did this when computing
% the houlry cloud fraction map.
% prescribe NaN if there is only 1 cloud fraction map on a given day.
% (e.g., Feb 16);

BndBox.cen = blobsIn.GeoLocs;
BndBox.Size = blobsIn.BoundingBoxSize;  
EqvDiam_dgr = blobsIn.stats_selected.EqvDiam .* blobsIn.L4_xres/111E3;  % in degrees.
search_RadRatio = search_box.sRadiusRatio;
cutout_extend = search_RadRatio;
for ib = 1:length(search_box.BndBox_dw)
    % get data in the search box

    %% how to do this when the search area is a circle?
    % rewrite the code below, so that I can choose how to select the region
    % of interest. 
    % use circle as the search area instead.
     % use the center of the blob and the distance for searching.
     
     xb = BndBox.cen(ib,1); yb = BndBox.cen(ib,2);
     if strcmp(masktype, 'circle') % do the following
         
         L = max(BndBox.Size(ib,:))*search_RadRatio;
         
         dist = sqrt((LON_cldfrac - xb).^2 + (LAT_cldfrac-yb).^2);
         smask = dist<=L*1.5;
         search_lon = LON_cldfrac(smask); search_lat = LAT_cldfrac(smask);
         
         search_lon = LON_cldfrac(smask);
         search_lat = LAT_cldfrac(smask);
         
         % reference cloud fraction: 2month large scale background mean;
         % small scale removed.
         cldInfo.refCF_LP2monthAve = refCF(smask);
         
         % reference cloud fraction: daily large scale background mean;
         cldInfo.refCF_LPdaily = cldfracbg_today(smask);
         
         % extract daily averaged cloud fraction at selected location:
         cldInfo.CF = cldfrac_today(smask);
         
         % absolute anomaly;
         cldInfo.CF_anom  = cldfrac_highpassed_today(smask);
                  
         % conundrum: how to define anomaly?  % does this make any sense here??
         cldInfo.CF_relanom_daily = cldInfo.CF_anom./cldInfo.refCF_LPdaily;
         cldInfo.CF_relanom_2mon = cldInfo.CF_anom./cldInfo.refCF_LP2monthAve;   % the anomaly here contains different information:
        
         
         
     elseif strcmp(masktype, 'square')
         
         L = EqvDiam_dgr(ib)/2*search_RadRatio;
         
         %% under construction:
         % use the same way to cut out cloudiness data as the SST data for
         % comparison.
         LX = BndBox.Size(ib,1);
         LY = BndBox.Size(ib,2);
         
         ExtRad = EqvDiam_dgr(ib)/2*cutout_extend ;
         
         lon0 = xb - LX/2; %- cutout_extend;   % index of the end of the box.
         lonN = xb + LX/2; %+ cutout_extend;
         
         lat0 = yb - LY/2; %- cutout_extend;
         latN = yb + LY/2; %+ cutout_extend;
         
         % create mask:
         cutout_mask = false(size(cldfrac_today));
         cutout_mask((LON_cldfrac>(lon0-ExtRad))&(LON_cldfrac<(lonN+ExtRad)))= true;
         cutout_mask((LAT_cldfrac<(lat0-ExtRad))|(LAT_cldfrac>(latN+ExtRad))) = false;
         
         latslice = LAT_cldfrac(:,1);
         nrow = length(find((latslice>=lat0-ExtRad)&(latslice<=latN+ExtRad)));
         
         % extract the pixels from ref_CF:
         %refCF_cutouts = refCF_today(cutout_mask);
         refCF_cutouts = refCF(cutout_mask);
         refCF_cutouts2D = reshape(refCF_cutouts, nrow, []);
         cldInfo.refCF_LP2monthAve = refCF_cutouts2D;
         
         dailyCFbg_cutouts = cldfracbg_today(cutout_mask);
         dailyCFbg_cutouts2D = reshape(dailyCFbg_cutouts, nrow, []);
         cldInfo.refCF_LPdaily = dailyCFbg_cutouts2D;
         
         
         % extract daily averaged cloud fraction at selected location:
         cldfrac_cutouts = cldfrac_today(cutout_mask);
         cldInfo.CF = reshape(cldfrac_cutouts, nrow, []);
         
         cldfrac_highpassed_cutouts = cldfrac_highpassed_today(cutout_mask);
         %cldInfo.CF_anom = (reshape(cldfrac_highpassed_cutouts, nrow, []))./refCF_cutouts2D;
         
         % absolute anomaly;
         cldInfo.CF_anom = reshape(cldfrac_highpassed_cutouts, nrow, []);
         
         %cldInfo.CF = cldfrac_today(smask);
         
         % conundrum: how to define anomaly?  % does this make any sense here??
         cldInfo.CF_relanom_daily = cldInfo.CF_anom./dailyCFbg_cutouts2D;
         cldInfo.CF_relanom_2mon = cldInfo.CF_anom./refCF_cutouts2D;   % the anomaly here contains different information:
         
         
         
         % if strcmp(masktype, 'square')
         search_lonvec = LON_cldfrac(cutout_mask); search_lon = reshape(search_lonvec, nrow, []);
         search_latvec = LAT_cldfrac(cutout_mask); search_lat = reshape(search_latvec, nrow, []);
         
         
     end
     

     
    
    if nhrly_maps<12       % consider these cloudiness frequency and anomaly value as invalid b/c the cloud mask sample size is too small.
        cldInfo.CF = NaN;
        cldInfo.CF_anom = NaN;
        cldInfo.CF_relanom_daily = NaN;
        cldInfo.CF_relanom_2mon = NaN;
        cldInfo.refCF_LP2monthAve = NaN;
        cldInfo.refCF_LPdaily = NaN;
        
    end
        

        
%     if length(time_all)>=6
%         cldfreq = compute_cloudfreq(SST_cloudmask);
%     else
    
    
    % visualize:
%     if ~isnan(cldInfo.CF) & ~isempty(cldInfo.CF)
%         BlobMask = blobsIn.blob_image{ib};
%         BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);
%         %CF_all = compute_cloudfreq(SST_all);
%         
%         figure(7);
%         if ib ==1
%             %pcolor(GOES_LON, GOES_LAT, CF_all); shading flat;
%             pcolor(LON_cldfrac, LAT_cldfrac, cldfrac_highpassed_today); shading flat; hold on
%             %scatter(search_lon(:), search_lat(:), 30, cldInfo.CF(:),'filled','MakerEdgeColor','k');
%         end
%         hold on;
%         % plot the warm blob element contours on here:
%         contour(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask),'-m','linewidth',1.2);
%         %hold off
%         colormap(jet);
%         colorbar;%wav
%         %pause(0.1)
%     end
%     
    % save data: 
    cloudInfo(ib).time = t;
    cloudInfo(ib).num_of_hrly_maps = nhrly_maps;
    cloudInfo(ib).cloudfrac = cldInfo.CF;
    cloudInfo(ib).cloudfrac_highfreq = cldInfo.CF_anom;
    cloudInfo(ib).cloudfrac_relanom_d = cldInfo.CF_relanom_daily;
    cloudInfo(ib).cloudfrac_relanom_2mon = cldInfo.CF_relanom_2mon;
   
    cloudInfo(ib).refCF_LP2monmean = cldInfo.refCF_LP2monthAve;
    cloudInfo(ib).refCF_LPdaily = cldInfo.refCF_LPdaily;
    cloudInfo(ib).lon = search_lon;
    cloudInfo(ib).lat = search_lat;
    
    
end


return