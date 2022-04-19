function cloudInfo = compute_cloudiness_over_features_from_L3C_cloudmasks(t, search_box, blobsIn, refCF_matfile, refstate)
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


% load in the reference cloudiness:
CF_refdata = load(['./cloudiness_data_store' filesep refCF_matfile]);
if strcmpi(refstate, 'JanFebMean')
    ref_CF = CF_refdata.CF_JanFeb;
elseif strcmpi(refstate, '3d_movmean')
    ref_CF = CF_refdata.CF_movmean_2monthsAve.wdsz_3d;
elseif strcmpi(refstate, '5d_movmean')
    ref_CF = CF_refdata.CF_movmean_2monthsAve.wdsz_5d;
else
    disp('invalid refstate, please enter either *JanFebMean*, *3d_movmean*, or *5d_movmean*. For more options, please generate the reference state first!')
    return
end


L3C_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/g16';
% load in L3C daily dataset;
tday = floor(t);
dateID = datestr(tday, 'mmmdd');
DOY = tday - datenum(2019, 12, 31);

% daily cloud mask.
matFN = ['GOES_SST_' num2str(DOY,'%3.3i') '_' dateID '_lon62W_TO_42W_lat5N_TO_25N.mat'];
disp(['loading' dateID])
load([L3C_datadir filesep matFN]);

[GOES_LON, GOES_LAT] = meshgrid(GOES_ATOMIC.lon(:,1), GOES_ATOMIC.lat(:,1));

% To-Consider: add in an option to select the level of clear-sky or cloud mask
QC_level = GOES_ATOMIC.quality_level;
% if level is selected, then this a selection criterion will be applied to
% this matrix.

SST_all = permute(GOES_ATOMIC.sea_surface_temperature, [2,1,3]);    % this is already the best_quality
time_all = GOES_ATOMIC.time_num;

% check cloud mask and throw away maps that contains all cloudy pixels.
% (I am considering it bad data)
%% QC data:
% throw away maps that is 100% garbage; How many of them in a day?)
nt0 = size(SST_all,3); nt = nt0;
i=1;
while i<=nt
    SST_tmp = SST_all(:,:,i);
    idxs = find(isnan(SST_tmp));
    if length(idxs) == numel(SST_tmp) 
        SST_all(:,:,i) = [];
        time_all(i) = [];
        
        nt = size(SST_all,3);
    else
        i = i+1;
    end
    
end

BndBox.cen = blobsIn.GeoLocs;
BndBox.Size = blobsIn.BoundingBoxSize;  
EqvDiam_dgr = blobsIn.stats_selected.EqvDiam .* blobsIn.L4_xres/111E3;  % in degrees.
search_RadRatio = search_box.sRadiusRatio;

for ib = 1:length(search_box.BndBox_dw)
    % get data in the search box
    if 1==0
        vertices = search_box.SBox_dw_vertices{ib};
        search_square_vx = vertices(1,:);
        search_square_vy = vertices(2,:);
        
        [in, on] = inpolygon(GOES_LON, GOES_LAT, search_square_vx, search_square_vy);
        smask = in | on;
    end
    
    %% how to do this when the search area is a circle?
     % use the center of the blob and the distance for searching.
    % L = max(BndBox.Size(ib,:))*search_RadRatio;
     L = EqvDiam_dgr(ib)/2*search_RadRatio;
     xb = BndBox.cen(ib,1); yb = BndBox.cen(ib,2);
     dist = sqrt((GOES_LON - xb).^2 + (GOES_LAT-yb).^2);
     smask = dist<=L*1.5;

    % extract the pixels from ref_CF:
    refCF_cutouts = ref_CF(smask);        

    % extract cloud mask at selected location for the entire day:
    SST_cloudmask= [];
    for it =1:length(time_all)
        SST_tmp = SST_all(:,:,it);
        SST_cloudmask(:,it) = SST_tmp(smask);
    end
    search_lon = GOES_LON(smask); search_lat = GOES_LAT(smask);
    
    % compute cloudiness (6:3:24) for different period and see how it
    % evolves. ?
    
    % here, I can call a different function to compute the cloudiness
    % anomalies. I will need a reference state though.
    % the good thing about this function is that the SST maps are
    % preserved.
    
       
    cldInfo = compute_cloudfreq_and_anom(SST_cloudmask, refCF_cutouts);
    
%     if length(time_all)>=6
%         cldfreq = compute_cloudfreq(SST_cloudmask);
%     else
    if length(time_all)<12       % consider these cloudiness frequency and anomaly value as invalid b/c the cloud mask sample size is too small.
        cldInfo.CF = NaN;
        cldInfo.CF_anom = NaN;
    end
    
    % visualize:
%     if ~isnan(cldInfo.CF)
%         BlobMask = blobsIn.blob_image{ib};
%         BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);
%         CF_all = compute_cloudfreq(SST_all);
%         
%         figure(7);
%         if ib ==1
%             pcolor(GOES_LON, GOES_LAT, CF_all); shading flat;
%         end
%         hold on;
%         % plot the warm blob element contours on here:
%         contour(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask),'-w','linewidth',1.2);
%         %hold off
%         colormap(jet);
%         colorbar;
%         %pause(0.1)
%     end
    
    % save data: 
    cloudInfo(ib).SST_cloudmask = SST_cloudmask;    % since the cloudmask is saved, I can then use it to find the cloudiness on a coarse grid.
    cloudInfo(ib).time = time_all;
    cloudInfo(ib).cloudfreq = cldInfo.CF;
    cloudInfo(ib).cloudfreq_anom = cldInfo.CF_anom;
    cloudInfo(ib).lon = search_lon;
    cloudInfo(ib).lat = search_lat;
    
    
end


return