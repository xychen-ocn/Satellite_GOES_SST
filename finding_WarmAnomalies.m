%%%--------------------------------------------------------------------------%%%
% Purpose: This script will be used to generate warm blob statistics using different
%          estimate of the temporal mean SST. Results will be saved to
%          blob_data under different folder names.
%
% Code Structure: 
%  - specify parameters to define the calculation, and detection of warm anomalies
%  - gather all the available g16 L3C data
%  - loop through the each daily dataset, determine if there is a large
%    scale mean field associated with it, if so, get the large scale mean
%    field out for use.
%  - compute the anomaly field of SST and send it to the detection
%    algorithm.
%  - if it is a RHB day of interest, then interpolate g16 data on the RHB
%    trajectory (and compute the SST anomaly from RHB for anomaly comparison)
%  - visualize the detected warm blobs (with RHB information imbedded if
%    it is a RHB day)
%  - save all warm blob statistics and warm blobs sampled by RHB along the
%    segments of interest.
% 
% Date: drafted on Jan 30, 2022 (XYC)
%
%%%--------------------------------------------------------------------------%%%

clear all; clc; close all;
addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');

global LON LAT
%% 0. load in all the needed datasets
load('RHB_wavelet_coherent_segment_Info.mat', 'RHB_segtime');
RHB_segdays = unique(cellstr(datestr(RHB_segtime,'mmmdd')));                                 % include days within the segment time.

% load a version of the RHB data that has the full RHB trajectories:
RHB_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
% load dates denote RHB SST variance:
load([RHB_datadir filesep 'rhb_daily_grouped_10min_data_0909latest.mat'],'rhbdates');

% load RHB full dataset:
load([RHB_datadir filesep 'EUREC4A_ATOMIC_RHB_AllinOne_v1.3.mat']);
RHB_full.time = RHB.time;
RHB_full.lon = RHB.lon; RHB_full.lat = RHB.lat;

clear RHB  % to reduce the amount of data loaded for speed and safety.

strong_SSTvar_days = cellstr(datestr(rhbdates.strong_SSTvar,'mmmdd'));
% broaden up the RHB_days to the 8 days with slightly more significant SST
% variance along with days with qualified segments examined.
RHB_days = unique(vertcat(RHB_segdays, strong_SSTvar_days));

% RHB_days = {'Jan01'};
% RHB_segdays = {'Jan01'};


% parameters for object (warm blobs) detection:
SSTgap_thres = 0.5;         % when cloudy pixel is less than 50% of the scene, fill the SST map;
dT_thres=0.1;               % temperature anomaly threshould.
cldcv_thres = 0.1;          % threshold to constrain whether or not a warm blob identified is under strong influence of cloud overhead (a.k.a., spatial interpolation.).
                            % threshold allow warm blobs that has very little cloudy pixel influence to
                            % be kept as the good warm blobs. 
minArea=25;                 % 5x5 pixels: 10km (1pixel is 2km across)
dT_thres_arry=[0:0.1:0.8];        % use different threshold to get all the blobs out.

meanSST_Options= {'g16','gpblend','ostia'};
AveWindows = [3, 5, 7, 15, 30];          % units: day
SST_AveMethod = {'moving','moving','moving','fixed','fixed'};              % 'moving': moving averaged; 'fixed':average during a fixed window length

search_region = [-59, -45; 8 22];    % avoid fft complication by land during 
keep_region = [-58, -48; 8 18];

LPF_scale_km = 600;                   % low-pass filtered scale (above which, the temporal mean SST will be subtacted from the total field.)

% plotting parameters:
RHBcolor = [0.8500 0.3250 0.0980];

%% 1. gather all the available g16 L3C data:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
files = dir([datadir filesep 'g16' filesep 'GOES_SST*.mat']);
filenames = {files.name};
% for i = 1:length(filenames)
%     g16FN = filenames{i};
%     load([datadir filesep 'g16' filesep g16FN]);
%     L3C_data(i) = GOES_ATOMIC;
% end

%% 2. ready to loop through different test parameters:
LPF_flag = true;
if LPF_flag
    tag2 = ['LPF' num2str(LPF_scale_km) 'km'];
else
    tag2 = ['no_spatial_filtering'];
end

tag3 = []; %'no_cloudmask';
for k = 1  %:2 %:length(meanSST_Options)
    product = meanSST_Options{k};
    clear SST_anom_along_RHBtraj
    
    for j = 2 %1:3  %length(AveWindows)
        wz = AveWindows(j);
        tag1 = [num2str(wz) 'd-' SST_AveMethod{j} 'Mean'];
        svtag =['wndsz_', num2str(wz) 'd'];
        
        SST_mean = get_large_scale_temporal_mean(product, wz, LPF_scale_km, search_region);
                
        %% build up folders to save data and figures:
        if ~isempty(tag3)
            casename = strjoin({product, tag1, tag2, tag3 },'_');
        else
            casename = strjoin({product, tag1, tag2 },'_');
        end

        datasvdir = [datadir filesep 'blob_data' filesep casename];
        if exist(datasvdir, 'dir')==0
            mkdir(datasvdir)
        end
        
        figsvdir = [datasvdir filesep 'figs'];
        if exist(figsvdir, 'dir')==0
            mkdir(figsvdir)
        end
        
        %% ready to loop through the daily dataset to find the warm anomalies:
        cnt_RHB = 0; Nsample = zeros(1,length(filenames));
        for i = 9; %1:length(filenames)
            g16FN = filenames{i};
            load([datadir filesep 'g16' filesep g16FN]);
            
            times_here = GOES_ATOMIC.time_num;
            tday = floor(times_here(1));
            dateID = datestr(tday, 'mmmdd');
            isRHB_day = ismember(dateID, RHB_days);
            
            blob_statFN = [dateID '_individual_warm_blobs_mapSelCCthres' num2str(SSTgap_thres)  ...
                '_blobSelCCthres' num2str(cldcv_thres) '.mat'];
            
            if true %~exist([datasvdir filesep blob_statFN])
                
                [SSTLS_available,meanSST,meanSST_LS] = check_availability_of_LS_meanSST(SST_mean, times_here);   % if available, the low-pass filtered temporal mean SST will be returned as well.
                %% do not apply spatial filtering: only do the temporal filtering at this point:
                
                % detect warm blobs if there is a large scale temporal mean field available:
                if SSTLS_available 
                    
                    disp(['--> working on ' dateID]);
                    
                    
                    
                    %% take a subset to search for blobs:
                    % take subset to exclude islands...
                    % longitude mask below actually applies to all the time instance:
                    lonmask = GOES_ATOMIC.lon(:,1)>=search_region(1,1) & GOES_ATOMIC.lon(:,1)<=search_region(1,2);
                    latmask = GOES_ATOMIC.lat(:,1)>=search_region(2,1) & GOES_ATOMIC.lat(:,1)<=search_region(2,2);
                    datasub.lon = double(GOES_ATOMIC.lon(lonmask,1));
                    datasub.lat = double(GOES_ATOMIC.lat(latmask,1));
                    
                    [LON, LAT]= meshgrid(datasub.lon, datasub.lat);
                    
                    dx = mean(mean(diff(datasub.lon)).*111.*cosd(datasub.lat));
                    dy = mean(diff(datasub.lat))*111;
                    spatialres = mean([dx,dy]);
                    
                    %% compute SST anomaly and filled the data gap to prep for blob detection.
                    SST_raw = permute(GOES_ATOMIC.sea_surface_temperature(lonmask,latmask,:), [2,1,3]);
                    [SST_filled, filled_flag] = fill_cloudy_pixels(SST_raw, LON, LAT, times_here, SSTgap_thres);
                    
                    if any(size(LON)~=size(meanSST))
                        % interpolate meanSST to the size of LON, LAT;
                        % (coarse --> finer) --> some spatial smoothing
                        % happended already.
                        meanSST = interp2(SST_mean.lon, SST_mean.lat, meanSST, LON, LAT);
                        meanSST_LS = interp2(SST_mean.lon, SST_mean.lat, meanSST_LS, LON, LAT);
                        
                    end
                        
                    
                    if LPF_flag
                        SST_anom_filled = SST_filled - meanSST_LS;
                        SST_anom_raw = SST_filled - meanSST_LS;
                        
                    else % only temporal filtered:
                        SST_anom_filled = SST_filled - meanSST;
                        SST_anom_raw = SST_filled - meanSST;               % 
                    end
                    
                    %% find anomaly along RHB trajectory
                    if isRHB_day
                        cnt_RHB = cnt_RHB + 1;
                        % for later analysis use.
                        sa_tmp = get_satellite_SST_anom_along_RHBtraj(SST_anom_raw, LON, LAT, times_here, RHB_full);
                        SST_anom_along_RHBtraj(cnt_RHB).values.(svtag) = sa_tmp.values;
                        SST_anom_along_RHBtraj(cnt_RHB).lon = sa_tmp.lon;
                        SST_anom_along_RHBtraj(cnt_RHB).lat = sa_tmp.lat;
                        SST_anom_along_RHBtraj(cnt_RHB).time = sa_tmp.time;
                    end
                    
                    
                    %% loop through all the available time:
                    % intialize variables:
           
                    cnt = 0; cnt_RHBseg=0;  VX_save = []; VY_save=[];
                    clear images indiv_blobs_collection RHB_sampled_blobs sampledblob_bank
                    
                    for it = 1:length(times_here)
                        %% detect warm blobs at each time intance:
                        SST_map_now = squeeze(SST_raw(:,:,it));
                        cloud_flag = isnan(SST_map_now);
                        
                        % get RHB sampled SST and trajectory:
                        t_now = times_here(it);
                        cond_time = (t_now-RHB_segtime(:,1)).*(t_now-RHB_segtime(:,2));
                        isRHB_segtime = any(cond_time<=0);
                        
                        disp(datestr(t_now));
                        %% check if there is a filled map available:
                        if filled_flag(it)
                            
                            
                            SST_anom_now = squeeze(SST_anom_filled(:,:,it));
                            
                            clear blob_all
                            thres_cnt=0;
                            for ithres=1:length(dT_thres_arry)
                                dT_thres = dT_thres_arry(ithres);
                                
                                % key function to detect warmer objects according to 8-connectivity and some thresholds.
                                if strcmp(tag3,'no_cloudmask')
                                    fake_cloud_flag = isnan(SST_anom_now); 
                                    cloud_flag = fake_cloud_flag;
                                end
                                
                                if mod(it,20)==0 && ithres==2
                                    checkflag=true;
                                else
                                    checkflag=false;
                                end
                                
                                [blob_info, exitflag] = find_SST_blobs(SST_anom_now, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
                                    'LON_grid', LON, 'LAT_grid', LAT,'checkflag',checkflag);
                                
                                %Nblobs = length(blob_image);
                                
                                if exitflag>0
                                    thres_cnt = thres_cnt+1;
                                    blob_all(thres_cnt) = blob_info;
                                end
                                %pause(0.1)
                                
                            end
                            
                            if thres_cnt > 0
                                cnt = cnt+1;
                                % call a function to remove the same captures of a SST blob:
                                debug=false;
                                [individual_blobs, Nblobs] = remove_doppelgangers(blob_all,debug);
                                individual_blobs.time = t_now;
                                
                                SST_filled_now = squeeze(SST_filled(:,:,it));
                                % should I feed in the SST anomaly instead of the
                                % filled SST map?
                                if LPF_flag
                                    indiv_blobs_collection(cnt) = retrieve_SST_and_relatedInfo_for_each_warm_blob(individual_blobs, ...
                                        LON, LAT, SST_filled_now, meanSST_LS);
                                else
                                    indiv_blobs_collection(cnt) = retrieve_SST_and_relatedInfo_for_each_warm_blob(individual_blobs, ...
                                        LON, LAT, SST_filled_now, meanSST);
                                end
                                    
                                % from this step I get SST_cutouts, SST_anom_cutouts, as well
                                % as the maximum SST anomaly at this location, and the
                                % background SST associated with this individual blob.
                                
                                
                                %% make figure to visualize the individual blobs caputure at this time instance (it):
                                if isRHB_day
                                
                                
                                    
                                    %% find blob sampled by RHB
                                    % the RHB data will be fetched by the function below (from
                                    % the object oriented scripts.)
                                    % the following function needs to be improved...to
                                    % allow more field to store in results.
                                    if isRHB_segtime
                                        
                                        hfig = visualize_selected_blobs(LON, LAT, SST_anom_now, indiv_blobs_collection(cnt), cloud_flag);
                                       %hfig = visualize_selected_blobs(LON, LAT, SST_anom_now+meanSST_LS-273.15, indiv_blobs_collection(cnt), cloud_flag, 'rainbow', 64);
                                       title(datestr(t_now),'fontsize',14);
                                       
                                       
                                        [results, RHB_data, found_blobs ]= search_for_blobs_sampled_by_platform('RHB', indiv_blobs_collection(cnt));
                                        if found_blobs
                                            cnt_RHBseg = cnt_RHBseg+1;
                                            RHB_sampled_blobs(cnt_RHBseg) = results;
                                        end
                                        
                                        tmask = (RHB_data.time <= t_now);
                                        [~,tnowID] = min(abs(RHB_data.time -t_now));
                                        
                                        % highlight the RHB (or other platform) sampled blobs:
                                        %%%%%% reserve code block %%%%%%%%%
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        hold on;
                                        plot(RHB_full.lon, RHB_full.lat, ':k','linewidth',1.2);  % full trajectory
                                        plot(RHB_data.lon, RHB_data.lat, '--k','linewidth',1.1);   % on the segment of interests.
                                        plot(RHB_data.lon(tmask), RHB_data.lat(tmask),'color', RHBcolor,'linewidth',1.2);
                                        hold on;
                                        plot(RHB_data.lon(tnowID), RHB_data.lat(tnowID),'color', RHBcolor, 'marker','p','MarkerFaceColor', RHBcolor, 'MarkerEdgeColor','k', 'MarkerSize',12);      % plot RHB trajectory
                                        if found_blobs
                                            plot(results.GeoLocs(:,1), results.GeoLocs(:,2), '*r','MarkerSize',10,'linewidth',1.2);
                                            % plot bounding box:
                                            %                                 [VX, VY] = get_boundingbox_vertices(results);
                                            %                                 plot(VX, VY, '-m','linewidth',1.5);
                                            % plot contour instead:
                                            for ib = 1:length(results.blobImage)
                                                contour(results.blobImageCoord(ib).lon, results.blobImageCoord(ib).lat, results.blobImage{ib},'-r','linewidth',1.);
                                            end
                                            
                                            sampledblob_bank(cnt_RHBseg) = results;
                                            
                                        end
                                        
                                        if cnt_RHBseg>1
                                            for ii = 1:cnt_RHBseg-1
                                                % plot the previously detected blobs again:
                                                plot(RHB_sampled_blobs(ii).GeoLocs(:,1), RHB_sampled_blobs(ii).GeoLocs(:,2),'*m');
                                                
                                            end
                                            prev_blobs_images = [sampledblob_bank.blobImage];
                                            prev_blobs_coord = [sampledblob_bank.blobImageCoord];
                                            
                                            for ib = 1:length(prev_blobs_coord)
                                                contour(prev_blobs_coord(ib).lon, prev_blobs_coord(ib).lat, prev_blobs_images{ib},':m','linewidth',0.8);
                                            end
                                        end
                                        
                                    else
                                        % not segment time, but still can show where RHB is sampling:
%                                         hold on;
%                                         plot(RHB_full.lon, RHB_full.lat, ':k','linewidth',1.2);  % full trajectory
%                                         
%                                         tmask = (RHB_full.time <= t_now)&(RHB_full.time>=times_here(1));
%                                         [~,tnowID] = min(abs(RHB_full.time -t_now));
%                                         plot(RHB_full.lon(tmask), RHB_full.lat(tmask),'color', RHBcolor,'linewidth',1.2);
%                                         hold on;
%                                         plot(RHB_full.lon(tnowID), RHB_full.lat(tnowID),'color', RHBcolor, 'marker','p','MarkerFaceColor', RHBcolor, 'MarkerEdgeColor','k', 'MarkerSize',12);      % plot RHB trajectory
%                                         
                                    end
                                    
                                    %caxis([-1, 1]);
                                    axis('equal');
                                    xlim([-59, -48]);
                                    ylim([8, 18]);
                                    
                                    
                                    %% save the daily results into an animation.
                                    % save frame for gif:
                                    frame = getframe(hfig);
                                    images{cnt} = frame2im(frame);
                                    
                                end
                               
                                
                            end                         % end thres_cnt >0
                            
                        end                             % add filled flag ==1
                        
                    end                                 % finished processing daily records.
                    
                    %% save data & animation if there is any blobs found on the current day:
                    if cnt > 0
                        % collect blob statistics in a better format:
%                         indiv_blobs_stats = get_blob_characters(indiv_blobs_collection, spatialres);
%                         
%                         blob_statFN = [dateID '_individual_warm_blobs_mapSelCCthres' num2str(SSTgap_thres)  ...
%                             '_blobSelCCthres' num2str(cldcv_thres) '.mat'];
%                         
%                         disp(['Saving detected blob info to ' datasvdir '...']);
%                         save([datasvdir filesep blob_statFN], 'indiv_blobs_collection', 'indiv_blobs_stats');
%                         
                        if isRHB_day
                            % get statistics for the RHB sampled blobs and averaged RHB sampled blobs.
%                             if cnt_RHBseg>1
%                                 RHB_blobs_stats = get_blob_characters(RHB_sampled_blobs, spatialres);
%                                 RHBblob_statFN = [dateID '_RHB_sampled_blobs.mat'];
%                                 
%                                 disp(['Saving blobs sampled by RHB info to ' datasvdir '...']);
%                                 save([datasvdir filesep RHBblob_statFN],'RHB_sampled_blobs','RHB_blobs_stats');
%                             end
                            
                        
                        %
                        
                        %% save animation:
                      % gifname = [dateID '_individual_warm_blobs_identified_' date '_with_SSTbg.gif'];

                        gifname = [dateID '_individual_warm_blobs_identified_' date '_redblue_colormap.gif'];
                        gif_abspath = [figsvdir filesep gifname];
                        write_frames_into_video(images, gif_abspath, 'gif');
                        end
                        
                        Nsamples(i) =  cnt;
                        
                    end
                    
                    
                    
           
                    
                end                               % operatable data:
            end                               % data already existed.
            
            
        end                                   % finish looping through daily data
        
        
        
    end                                        % loop through different cases
    
    % save the anomaly data computed during all the RHB days with this dataproduct:
%     anom_FN = [product 'SST_anomaly_extracted_along_RHBtraj.mat'];
%     save([datasvdir filesep anom_FN], 'SST_anom_along_RHBtraj','RHB_days');

end
