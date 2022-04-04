%% Purpose: This script is used to find warm patches from L4 daily SST products.
%% these warm patches are the small scale spatial SST anomalies. 
% 1. Code Structure:
%    - find warm anomalies
%    - get blob statistics
%    - get averaged ERA5 wind direction within 2x bounding box?
%    - get features-bounding box centric coordinate with its x-axis aligned
%    with the wind direction. 
%    - compute daily cloudiness above the feature in the coordinate defined
%    above (using L3C SST cloud mask --> need a utility function for this.)
% 
% 2. References:
%  the idea comes from Park et al. (2006)
%  part of the codes from: finding_WarmAnomalies.m
%
% 3. Author:
% XChen (drafted at Feb 21, 2022);
% To-do: 
%  - sensitivity test for the SST threshold. (check how the size and number
%   and shape of the blobs are changed by the SST detection threshold)
%  - add background 2-month cloudiness as one of the output parameter. 

% To-fix:
% - 1. how the different types of day are selected
% (grouped_data_by_large_scale conditions)
% - 2. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

% load in L4 data
dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
productName = 'gpblend';    % is closer to g16 daily averaged SST than OSTIA.
L4matFN = [productName '_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N'];

load([dataroot filesep productName filesep L4matFN]);
L4data = data;  clear data

figsvdir = [dataroot filesep 'blob_data/L4test' filesep 'figs'];


% define parameters to find warm patches:
dT_thres=0.1;               % temperature anomaly threshould.
cldcv_thres = 0.1;          % threshold to constrain whether or not a warm blob identified is under strong influence of cloud overhead (a.k.a., spatial interpolation.).
minArea=25;                 % 5x5 pixels: 10km (1pixel is 2km across)
dT_thres_arry=[0.1:0.1:0.4];        % use different threshold to get all the blobs out.  we don't need this.

cutoff_scale_km = 600;

ERA5dataFN = 'era5_10mWindVectors_DJF2019-2022.nc';
search_RadRatio = 2;    % 2x the bounding box size. for wind speed

clear indiv_blobs_collection
for it = 1:length(L4data.time)
    %% a: find daily (averaged) SST warm patches:
    t_now = L4data.time_num(it);
    
    SST0 = L4data.analysed_sst(:,:,it)';
    lon = double(L4data.lon(:,it));
    lat = double(L4data.lat(:,it));
    
    L4res = mean([diff(lon)*111E3.*cosd(lat(1:end-1)), diff(lat)*111E3],'all');
    
    [LON, LAT] = meshgrid(lon, lat);
    
    if mod(it, 10)==0
        checkflag = true;
    else
        checkflag = false;
    end
    
    meanSST_LS =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);    % the size is different here.
    
    SST_anom_spatial = SST0 - meanSST_LS;
    fake_cloud_flag = isnan(SST_anom_spatial);
    cloud_flag = fake_cloud_flag;
    
    %     clear blob_all
    %     thres_cnt=0;
    %     for ithres=1:length(dT_thres_arry)
    %         dT_thres = dT_thres_arry(ithres);
    % dT_thres = 0.1
    
    % detect both cold and warm anomalies together. 
    [blob_info, exitflag] = find_SST_blobs(SST_anom_spatial, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
        'LON_grid', LON, 'LAT_grid', LAT,'checkflag', false);
    caxis([-0.8, 0.8]);
    colormap(redblue);
    
    %         if exitflag>0
    %             thres_cnt = thres_cnt+1;
    %             blob_all(thres_cnt) = blob_info;
    %         end
    
    
    
    %% b: get blob statistics
    %     if thres_cnt > 0
    %         debug = true;
    %         [individual_blobs, Nblobs] = remove_doppelgangers(blob_all, debug);
    %         individual_blobs.time = t_now;
    blob_info.time =t_now;
    blob_info.L4_xres = L4res;
    indiv_blobs_collection(it) = retrieve_SST_and_relatedInfo_for_each_warm_blob(blob_info, LON, LAT, SST0, meanSST_LS);
    
    % visualize it better and save as animations:
    if 1==0
        hfig = visualize_selected_blobs(LON, LAT, SST_anom_spatial, indiv_blobs_collection(it), cloud_flag);
        hold on;
        % plot the LS SST contours;
        [hc,c]= contour(LON, LAT, meanSST_LS-273.15, [25:0.5:28],'linestyle','--','linewidth',1.2,'color','b');
        clabel(hc,c,[25:0.5:28],'labelspacing',300,'color','b');
        
        title(datestr(t_now),'fontsize',14);
        
        xlim([-59, -48]);
        ylim([8, 18]);
        %axis('equal');
        caxis([-0.8, 0.8]);
        pause(0.1);
        hold off
        
        frame = getframe(hfig);
        images{it} = frame2im(frame);
    end
    
    if 1==1
        %% c: get averaged ERA5 wind direction within 2x bounding box area. (need a utility function here)
        %%    we could compute the surface wind divergence as well in the bounding box region. (need a utility function here)
        % to-do tomorrow: update ERA5 wind to include wind at the cloud
        % base and cloud top.
        [ave_windInfo, SearchArea_dw] = get_ERA5_windInfo_over_features(t_now, indiv_blobs_collection(it), search_RadRatio, ERA5dataFN);      % output: structure: dir, spd, u,v;
        
        
        %% d: get L3C cloud mask and compute cloudiness for the area within 2x bounding box. (need a utility function here)
        cloudiness = compute_cloudiness_over_features_from_L3C_cloudmasks(t_now, SearchArea_dw, indiv_blobs_collection(it));  % cloudiness is a structure, contain cloudiness over each bounding box.
        
        %% e: get wind aligned, feature centric, normalized coordinate to show the cloudiness map. (need a utility function here)
        % add option to approx feature by circle instead of by ellipse. 
        cloudiness_downwind_FCN{it} = map_cloudiness_to_newCoord(ave_windInfo,  cloudiness, indiv_blobs_collection(it));
        
        %indiv_blobs_collection(it).cloudiness_downwind_FCN = cloudiness_downwind_FCN;
        clear cloudiness
        %% f. compute the ERA5 wind convergence in the normalized coordinate. (can be done later)
        
        close all
    end
    %end
end

datasvdir = [dataroot filesep 'blob_data/L4test'];
if ~exist(datasvdir)
    mkdir(datasvdir)
end
save([datasvdir filesep 'L4_gpblend_daily_WarmSSTanom_blobs_and_cloudiness_larger_feature_domain.mat'],'indiv_blobs_collection','cloudiness_downwind_FCN','-v7.3');



gifname = ['JanFeb2020_L4_gpblend_individual_warm_blobs_identified_dTthres0.1_LPF600km_smallRegion.gif'];
gif_abspath = [figsvdir filesep gifname];
write_frames_into_video(images, gif_abspath, 'gif');


if 1==0
%% cold patches:
clear indiv_blobs_collection
for it = 1:length(L4data.time)
    %% a: find daily (averaged) SST warm patches:
    t_now = L4data.time_num(it);
    
    SST0 = L4data.analysed_sst(:,:,it)';
    lon = double(L4data.lon(:,it));
    lat = double(L4data.lat(:,it));
    
    L4res = mean([diff(lon)*111E3.*cosd(lat(1:end-1)), diff(lat)*111E3],'all');
    
    [LON, LAT] = meshgrid(lon, lat);
    
    if mod(it, 10)==0
        checkflag = true;
    else
        checkflag = false;
    end
    
    meanSST_LS =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);    % the size is different here.
    
    SST_anom_spatial = SST0 - meanSST_LS;
    fake_cloud_flag = isnan(SST_anom_spatial);
    cloud_flag = fake_cloud_flag;
    
    %     clear blob_all
    %     thres_cnt=0;
    %     for ithres=1:length(dT_thres_arry)
    %         dT_thres = dT_thres_arry(ithres);
    % dT_thres = 0.1
    
    [blob_info, exitflag] = find_SST_blobs(-SST_anom_spatial, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
        'LON_grid', LON, 'LAT_grid', LAT,'checkflag', checkflag);
    caxis([-0.8, 0.8]);
    colormap(redblue);
    
    %         if exitflag>0
    %             thres_cnt = thres_cnt+1;
    %             blob_all(thres_cnt) = blob_info;
    %         end
    
    
    
    %% b: get blob statistics
    %     if thres_cnt > 0
    %         debug = true;
    %         [individual_blobs, Nblobs] = remove_doppelgangers(blob_all, debug);
    %         individual_blobs.time = t_now;
    blob_info.time =t_now;
    blob_info.L4_xres = L4res;
    indiv_blobs_collection(it) = retrieve_SST_and_relatedInfo_for_each_warm_blob(blob_info, LON, LAT, SST0, meanSST_LS, 0.1,'cold');
    
    % visualize it better and save as animations:
    if 1==0
        hfig = visualize_selected_blobs(LON, LAT, SST_anom_spatial, indiv_blobs_collection(it), cloud_flag);
        hold on;
        % plot the LS SST contours;
        [hc,c]= contour(LON, LAT, meanSST_LS-273.15, [25:0.5:28],'linestyle','--','linewidth',1.2,'color','b');
        clabel(hc,c,[25:0.5:28],'labelspacing',300,'color','b');
        
        title(datestr(t_now),'fontsize',14);
        
        xlim([-59, -48]);
        ylim([8, 18]);
        %axis('equal');
        caxis([-0.8, 0.8]);
        pause(0.1);
        hold off
        
        frame = getframe(hfig);
        images{it} = frame2im(frame);
    end
    
    if 1==1
        %% c: get averaged ERA5 wind direction within 2x bounding box area. (need a utility function here)
        %%    we could compute the surface wind divergence as well in the bounding box region. (need a utility function here)
        [ave_windInfo, SearchArea_dw] = get_ERA5_windInfo_over_features(t_now, indiv_blobs_collection(it), search_RadRatio, ERA5dataFN);      % output: structure: dir, spd, u,v;
        
        
        %% d: get L3C cloud mask and compute cloudiness for the area within 2x bounding box. (need a utility function here)
        cloudiness = compute_cloudiness_over_features_from_L3C_cloudmasks(t_now, SearchArea_dw, indiv_blobs_collection(it));  % cloudiness is a structure, contain cloudiness over each bounding box.
        
        %% e: get wind aligned, feature centric, normalized coordinate to show the cloudiness map. (need a utility function here)
        % this
        cloudiness_downwind_FCN{it} = map_cloudiness_to_newCoord(ave_windInfo,  cloudiness, indiv_blobs_collection(it));
        
        %indiv_blobs_collection(it).cloudiness_downwind_FCN = cloudiness_downwind_FCN;
        clear cloudiness
        %% f. compute the ERA5 wind convergence in the normalized coordinate. (can be done later)
        
        close all
    end
    %end
end
datasvdir = [dataroot filesep 'blob_data/L4test'];
if ~exist(datasvdir)
    mkdir(datasvdir)
end
save([datasvdir filesep 'L4_gpblend_daily_ColdSSTanom_blobs_and_cloudiness_larger_feature_domain.mat'],'indiv_blobs_collection','cloudiness_downwind_FCN','-v7.3');
end
