% This script will call the "blob" script to do feature identification from
% the SST map. The features capture here will also compared to eddies
% identified from an eddy-tracking code pacakge.

% Notes: 
%  issues: had problem separating pixels from the missing area.  ( I think
%  I need to fill in the cloudy gap first.. Yes, --> resolved 01/11/2021)
%  find_SST_blobs is an adapation of the code for the cloud image.
%
%  To-DO/Test: define a mean that accounts for the large scale gradient at
%  different latitude bands. (attach to each warm blob (SST anomaly) a background SST value.)
%
%global LON LAT
% code design: 
clear all; clc; 
addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');

%% manually created the following data:
% RHB_segtime = [time_start_valid, time_end_valid];           
% save('RHB_wavelet_coherent_segment_Info.mat', 'RHB_segtime');
load('RHB_wavelet_coherent_segment_Info.mat', 'RHB_segtime');
RHB_days = unique(cellstr(datestr(RHB_segtime,'mmmdd')));                                 % include days within the segment time.

% load a version of the RHB data that has the full RHB trajectories:
RHB_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
load([RHB_datadir filesep 'rhb_daily_grouped_10min_data_0909latest.mat']);


% load in all images:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
files = dir([datadir filesep 'GOES_SST*.mat']);
filenames = {files.name};

%load([datadir filesep 'GOES_SST_Jan09.mat']);
datasvdir0 = [datadir filesep 'blob_data'];
method = 'planefit';
datasvdir = [datadir filesep 'blob_data/SSTdetrend_BY_' method];
if exist(datasvdir, 'dir')==0
    mkdir(datasvdir)
end

figsvdir = [datasvdir filesep 'figs'];
if exist(figsvdir, 'dir')==0
    mkdir(figsvdir)
end

% what should be the mean to subtract the daily SST data from? 
% Try daily mean first, which should remove diurnal cycle of SST variability.
% the remaining SST variability would come from other processes.

% specify parameters for the warm blob searching mission!
SSTgap_thres = 0.5;   % when cloudy pixel is less than 50% of the scene, fill the SST map;
dT_thres=0.1;         % temperature anomaly threshould.
cldcv_thres = 0.1;    % threshold to constrain whether or not a warm blob identified is under strong influence of cloud overhead (a.k.a., spatial interpolation.).
                      % threshold allow warm blobs that has very little cloudy pixel influence to
                      % be kept as the good warm blobs. 
minArea=25;   % 5x5 pixels: 10km (1pixel is 2km across)
% not as straight forward as I expected.. different dT threshold
% yields different answers. (those answers do not overlap)

dT_thres_arry=[0:0.1:0.8];        % note that 75 prctile ~0.25 on Jan9;  old threshould when the domain mean SST is used.[0.01, 0.1:0.1:0.8];

% call another script to compute the 2-month SST mean;
% in the following function, all the satellite during the 2-month period in
% 2020 will be read and used to compute the mean SST large scale gradient.
% two method can be used to determine the LSG: plane fit; Fourier spectrum
% subtraction.
% methodstr = 'spectrum';     % or 'planefit';  if spectrum, need to select a scale.
% [SST_LS_trend_spec,SST_maps_local_mean_monthly] = estimate_ATOMIC_period_mean_SST_LSG('method',methodstr, 'CutoffScale', 1000);
% [SST_LS_trend_fit,SST_maps_local_mean_monthly] = estimate_ATOMIC_period_mean_SST_LSG;
% 

DOIs = unique(cellstr(datestr(RHB_segtime,'mmmdd')));
%DOIs = {'Feb04','Feb07','Jan24'};
%warm_anom_prctiles=zeros(length(filenames), 4);
% for j = 1:length(DOIs)
fid = fopen('warning.txt','w');
 for j = 1:length(filenames)
    matFN = filenames{j};
    tmp = strsplit(matFN,'_');
    dateID = tmp{end}(1:5);

%        dateID = DOIs{j};
%        matFN = ['GOES_SST_' dateID '.mat'];
    
    
    blob_statFN = [dateID '_individual_warm_blobs_noLargeScaleSSTGrad_in_LC3_GOES-SST_mapSelCCthres' num2str(SSTgap_thres)  ...
        '_blobSelCCthres' num2str(cldcv_thres) '.mat'];
    
    %if ~exist([datasvdir filesep blob_statFN])
        
        load([datadir filesep matFN]);
        
        disp(['--> working on ' dateID]);
        
        %% check to see if this is a RHB day of interest:
        isRHB_day = ismember(dateID, RHB_days);
        
        % take subset to exclude islands...
        xstart = -59;
        % longitude mask below actually applies to all the time instance:
        lonmask = GOES_ATOMIC.lon(:,1)>xstart;
        datasub.lon = double(GOES_ATOMIC.lon(lonmask,1));
        datasub.lat = double(GOES_ATOMIC.lat(:,1));
        
        [LON, LAT]= meshgrid(datasub.lon, datasub.lat);
        
        dx = mean(mean(diff(datasub.lon)).*111.*cosd(datasub.lat));
        dy = mean(diff(datasub.lat))*111;
        spatialres = mean([dx,dy]);
        
        
        %% compute daily averaged SST and the mean SST to find warm blobs.
        %     all_SST = reshape(GOES_ATOMIC.sea_surface_temperature(lonmask,:,:), [],1);
        %     SST_daily_mean = mean(all_SST, 'omitnan');
        
        %% - interpolate the SST map to fill the cloud gap;
        SST_raw = permute(GOES_ATOMIC.sea_surface_temperature(lonmask,:,:), [2,1,3]);
        %SSTgap_thres = 0.5;  % when cloudy pixel is less than 50% of the scene, fill the SST map;
        [SST_filled, filled_flag] = fill_cloudy_pixels(SST_raw, LON, LAT,GOES_ATOMIC.time_num, SSTgap_thres);
        
        %if length(find(filled_flag==1))/length(filled_flag)>0.1
            
            %% compute regional daily averaged SST to remove large scale SST trend.
            % lat = 8~ 18N.
            % three different SST meridional zones, 2 SST zonal zones)
            % compute the anomaly in this way; return a matrix;
            SST_maps = permute(GOES_ATOMIC.sea_surface_temperature(lonmask,:,:), [2,1,3]);  % NLAT x NLON x NTIME
            % not readly to work...
            if strcmp(method, 'spectrum')
                [SST_LS_trend, SST_maps_local_mean] = estimate_LSG_SST(datasub.lon, datasub.lat, ...
                    SST_maps, 'method','spectrum','CutoffScale',5000);
            else
                [SST_LS_trend, SST_maps_local_mean] = estimate_LSG_SST(datasub.lon, datasub.lat, ...
                    SST_maps, 'method','planefit');
            end
            
            %% note: need to add capacity to test for sensitivity to the way the anomalies are computed.
            
            %     % simply compute the meridionally averaged SST:
            %     %SST_meridional_mean = mean(SST_maps_local_mean, 1, 'omitnan');
            %     valid = ~isnan(SST_maps_local_mean);
            %     C =planefit(LON(valid), LAT(valid), SST_maps_local_mean(valid));
            %     SST_trend_plane = C(1)*LON + C(2)*LAT + C(3);
            %
            %%% --- call a function to compute SST mean ---- %%%
            % compute_monthly_mean_SST_and_anomaly()
            
            figure(1);
            subplot(1,3,1)
            pcolor(LON, LAT, SST_maps_local_mean); shading flat; colorbar;
            axis('square');
            title({'Jan-Feb mean SST map'; '(excluded nan)'});
            
            subplot(1,3,2)
            pcolor(LON, LAT, SST_LS_trend); shading flat; colorbar;
            axis('square');
            title({'Jan-Feb large scale SST gradient' ;'(from valid data point)'});
            
            subplot(1,3,3)
            pcolor(LON, LAT, SST_maps(:,:,12) - SST_LS_trend); shading flat; colorbar
            axis('square');
            title({'SST anomaly';' at 12UTC'});
            
            pause(0.2);
            xc_savefig(gcf, figsvdir, [dateID '_warm_anom_obtain_steps_SSTdetrendBY_' method '.jpg'], [0 0 10 4]);
            %     %close gcf
            %
            %     SST_anom_zonally_detrend_all = SST_maps - repmat(SST_trend_plane, 1 , 1, nt);
            %
            %     figure(2); clf;
            %     histogram(SST_anom_zonally_detrend_all);
            %     xlabel('\Delta SST (K)');
            %     ylabel('count');
            %     grid on
            %     set(gca,'fontsize',12);
            %     title(dateID, 'fontsize',14);
            %     xc_savefig(gcf, figsvdir, [dateID '_warm_anom_distribution.jpg'], [0 0 10 8]);
            %     %close gcf
            %
            %     warm_anom = SST_anom_zonally_detrend_all(SST_anom_zonally_detrend_all>0);
            %     warm_anom_prctiles(j,:) = prctile(warm_anom,[25, 50, 75, 99]);
            
            %SST_LS_trend = SST_LS_trend_fit;
            cnt = 0; cnt_RHB=0; VX_save = []; VY_save=[];
            clear images indiv_blobs_collection RHB_sampled_blobs sampledblob_bank
            try
                for it =1:size(SST_maps,3)
                    SST_map_now = GOES_ATOMIC.sea_surface_temperature(lonmask,:,it)';
                    % mean_SST_hrly = mean(SST_map_now(:),'omitnan');
                    cloud_flag = isnan(SST_map_now);
                    
                    % get RHB sampled SST and trajectory:
                    t_SST = GOES_ATOMIC.time_num(it);
                    cond_time = (t_SST-RHB_segtime(:,1)).*(t_SST-RHB_segtime(:,2));
                    isRHB_segtime = any(cond_time<=0);
                    
                    disp(datestr(t_SST));
                    %% check if there is a filled map available:
                    if filled_flag(it)
                        cnt = cnt+1;
                        
                        % yes & carry out the search for warm blob:
                        SST_anom = squeeze(SST_filled(:,:,it)) - SST_LS_trend;
                        %SST_anom_daily = squeeze(SST_filled(:,:,it)) - SST_daily_mean;
                        %SST_anom_hrly = squeeze(SST_filled(:,:,it)) - mean_SST_hrly;
                        
                        %SST_anom_zonally_detrend = squeeze(SST_filled(:,:,it)) - SST_trend_plane;
                        
                        %SST_anom = SST_anom_zonally_detrend;
                        %SST_anom = SST_anom_daily;
                        
                        %for dT_thres = 0:0.1:0.6
                        clear blob_all
                        thres_cnt=0;
                        for ithres=1:length(dT_thres_arry)
                            dT_thres = dT_thres_arry(ithres);
                            
                            %
                            %[blob_coord, blob_image, SST_blobmask, blob_stats_tmp, exc_blob_stats] = ...
                            [blob_info, exitflag] = find_SST_blobs(SST_anom, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
                                'LON_grid', LON, 'LAT_grid', LAT);
                            
                            %Nblobs = length(blob_image);
                            
                            if exitflag>0
                                thres_cnt = thres_cnt+1;
                                blob_all(thres_cnt) = blob_info;
                                %                 else
                                %                     blob_all(ithres)=NaN;
                            end
                            pause(0.5)
                            
                        end
                        
                        if thres_cnt > 0
                            % call a function to remove the same captures of a SST blob:
                            debug=true;
                            [individual_blobs, Nblobs] = remove_doppelgangers(blob_all,debug);
                            individual_blobs.time = GOES_ATOMIC.time_num(it);
                            
                            SST_map_filled = squeeze(SST_filled(:,:,it));
                            % should I feed in the SST anomaly instead of the
                            % filled SST map?
                            indiv_blobs_collection(cnt) = retrieve_SST_and_relatedInfo_for_each_warm_blob(individual_blobs, ...
                                LON, LAT, SST_map_filled, SST_LS_trend);
                            % from this step I get SST_cutouts, SST_anom_cutouts, as well
                            % as the maximum SST anomaly at this location, and the
                            % background SST associated with this individual blob.
                            
                            
                            %% make figure showing the individual blobs caputure at this time instance (it):
                            
                            %   call a function to do this as well
                            hfig = visualize_selected_blobs(LON, LAT, SST_anom, indiv_blobs_collection(cnt), cloud_flag);
                            title(datestr(t_SST),'fontsize',14);
                            
                            %% save blobs that RHB sampled separately.
                            if isRHB_day %&& isRHB_segtime
                                % the RHB data will be fetched by the function below (from
                                % the object oriented scripts.)
                                % the following function needs to be improved...to
                                % allow more field to store in results.
                                if isRHB_segtime
                                    [results, RHB_data, found_blobs ]= search_for_blobs_sampled_by_platform('RHB', indiv_blobs_collection(cnt));
                                    if found_blobs
                                        cnt_RHB = cnt_RHB+1;
                                        RHB_sampled_blobs(cnt_RHB) = results;
                                    end
                                    
                                    tmask = (RHB_data.time <= GOES_ATOMIC.time_num(it));
                                    [~,tnowID] = min(abs(RHB_data.time -GOES_ATOMIC.time_num(it)));
                                    
                                    % highlight the RHB (or other platform) sampled blobs:
                                    %%%%%% reserve code block %%%%%%%%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    hold on;
                                    RHBcolor = [0.8500 0.3250 0.0980];
                                    plot(RHB.lon, RHB.lat, ':k','linewidth',1.2);  % full trajectory
                                    plot(RHB_data.lon, RHB_data.lat, '--k',1.1);   % on the segment of interests.
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
                                        
                                        sampledblob_bank(cnt_RHB) = results;
                                        
                                        % title({datestr(t_SST);[num2str(nRHBblobs) ' blobs sampled by RHB']},'fontsize',14);
                                    end
                                    
                                    if cnt_RHB>1
                                        for ii = 1:cnt_RHB-1
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
                                    hold on;
                                    plot(RHB.lon, RHB.lat, ':k','linewidth',1.2);  % full trajectory
                                    
                                    tmask = (RHB.time <= GOES_ATOMIC.time_num(it))&(RHB.time>=GOES_ATOMIC.time_num(1));
                                    [~,tnowID] = min(abs(RHB.time -GOES_ATOMIC.time_num(it)));
                                    plot(RHB.lon(tmask), RHB.lat(tmask),'color', RHBcolor,'linewidth',1.2);
                                    hold on;
                                    plot(RHB.lon(tnowID), RHB.lat(tnowID),'color', RHBcolor, 'marker','p','MarkerFaceColor', RHBcolor, 'MarkerEdgeColor','k', 'MarkerSize',12);      % plot RHB trajectory
                                    
                                end
                                
                                %caxis([-1, 1]);
                                xlim([-59, -48]);
                                ylim([8, 18]);
                                axis('equal');
                                
                                
                                
                            end
                            % save frame for gif:
                            frame = getframe(hfig);
                            images{cnt} = frame2im(frame);
                            
                        end
                        
                    end
                    
                    % pause
                end
                %end
                
                if cnt > 0
                    indiv_blobs_stats = get_blob_characters(indiv_blobs_collection, spatialres);
                    
                    % save any saved records
                    blob_statFN = [dateID '_individual_warm_blobs_in_LC3_GOES-SST_mapSelCCthres' num2str(SSTgap_thres)  ...
                        '_blobSelCCthres' num2str(cldcv_thres) '_SSTdetrendBY_' method '.mat'];
                    save([datasvdir filesep blob_statFN], 'indiv_blobs_collection', 'indiv_blobs_stats');
                    
                    if isRHB_day
                        % get statistics for the RHB sampled blobs and averaged RHB sampled blobs.
                        if cnt_RHB>1
                            RHB_blobs_stats = get_blob_characters(RHB_sampled_blobs, spatialres);
                            RHBblob_statFN = [dateID '_RHB_sampled_blobs_SSTdetrendBY_' method '.mat'];
                            save([datasvdir filesep RHBblob_statFN],'RHB_sampled_blobs','RHB_blobs_stats');
                        end
                    end
                    %
                    % save the identification process as as a gif:
                    %gifname = [dateID '_warm_blobs_identification_dTthres' num2str(dT_thres) '.gif'];
                    gifname = [dateID '_individual_warm_blobs_identified_detrendedSST_' method '.gif'];
                    gif_abspath = [figsvdir filesep gifname];
                    %gif_abspath0 = [datasvdir0 filesep gifname];
                    write_frames_into_video(images, gif_abspath, 'gif');
                    
                end
                
                %pause;
                %% average RHB sampled blobs on a given day:
                % get statistics for the RHB sampled blobs and averaged RHB sampled blobs.
                
                % strength (SST anomaly); Size (equiv. diameter); Major_Axis;
                % Minor_axis; Orientation;
                % this can be in another script.
                % these are blobs captured over the course of the segment.
                
                % now plot these statistics in a plot:
                % EqvDiam: diameter of a circle with the same area as the region.
                % (4*Area/pi)
                %             if isRHB_day && cnt_RHB>0
                %                 qtty_list = {'SST_anom','EqvDiam_km','AxisRatio', 'SST_bg','MajAxisLen','MinAxisLen'};
                %                 xlabels ={'SST anomaly (K)', 'Equivalent Diameter (km)', 'Axis Ratio', 'ave. background SST', 'Major Axis Length (km)','Minor Axis Length (km)'};
                %
                %                 figure(j+20);
                %                 set(gcf,'name',dateID);
                %                 for iq = 1:length(qtty_list)
                %                     qn = qtty_list{iq};
                %                     qval = indiv_blobs_stats.(qn);
                %                     RHB_qval = RHB_blobs_stats.(qn);
                %
                %                     subplot(2,3,iq);
                %                     if iq == 2 || iq==5 || iq ==6
                %                         binedges = [10:10:100];
                %                         histogram(qval, binedges);
                %                     else
                %                         histogram(qval);
                %                     end
                %                     hold on;
                %                     yrange = get(gca,'ylim');
                %                     plot(mean(qval)*ones(1,2),yrange,'-y','linewidth',1.2);
                %                     plot(median(qval)*ones(1,2), yrange,'-m','linewidth',1.2);
                %                     plot(prctile(qval,25)*ones(1,2), yrange, '--k');
                %                     plot(prctile(qval,75)*ones(1,2), yrange, '--k');
                %                     labelstr = strsplit(qn,'_');
                %                     xlabel(xlabels{iq});
                %
                %                     % add the mean val sampled by RHB:
                %                     plot(mean(RHB_qval)*ones(1,2), yrange, 'color', RHBcolor, 'linewidth',1.2);
                %
                %                     if iq == 2 || iq==5 || iq ==6
                %                         xlim([0 120]);
                %                     end
                %                     hold off
                %                     grid on
                %                     set(gca,'fontsize',14);
                %                 end
                %                 figname = [dateID '_blobs_statics_checkv0.jpg'];
                %                 xc_savefig(gcf, figsvdir, figname, [0 0 12 6]);
                %
                %             end
                %pause
            catch
                
                fprintf(fid,'%s\n',[dateID ' failed']);
                warning([dateID ' failed! needs to revisit.']);
            end
                
                    
 %end
    %end
    
end

close(fid);

%save([datasvdir filesep 'warm_anom_prctiles_25_50_75_99'],'warm_anom_prctiles');



%%%% Obsolete %%% ;
%                     figure(1); clf;
%                     subplot(1,2,1)
%                     pcolor(LON, LAT, SST_filled(:,:,it));
%                     shading flat;
%                     title('filled SST');
%                     colormap('parula')
%                     colorbar
%                     caxis([SST_daily_mean-1.5,SST_daily_mean+1.5]);
%                     axis('square')
%                     
%                     subplot(1,2,2)
%                     pcolor(LON, LAT, SST_anom_daily);
%                     shading flat;
%                     hold on
%                     contour(LON,LAT, SST_anom_daily, [-0.8:0.2:0.8],'k');
%                     title('SST anomaly');
%                     colormap('parula')
%                     colorbar
%                     caxis([-1.5,+1.5]);
%                     axis('square')
                    

