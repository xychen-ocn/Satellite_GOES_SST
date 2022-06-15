% This script will combine two previous script together to automate
% sensitivity testing (so that it becomes less of a mental drag for me).
% X.C. 6/6/2022
%
% what kind of output do I need to generate?
%  - a movie of the SST anomalies etc. 
%  - composite results
%  - individual case snapshots
% 
global dataroot casesvdir figsvdir datasvdir
global indiv_blobs_collection random_blobs
global downwind_FCN_regridded rb_downwind_FCN_regridded cutouts_physicalspace
global sample_random_blobs
% global cloudiness_downwind_FCN SSTcutout_downwind_FCN winddiv_downwind_FCN

% random blobs:
% global  cloudiness_downwind_FCN_rb 
% global SSTcutout_downwind_FCN_rb winddiv_downwind_FCN_rb


clear all; clc; close all;

%% 0. preparing directories:
dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
addpath(genpath([dataroot filesep 'util']));   
addpath(genpath([dataroot filesep 'probability_analysis']));

casesvdir = [dataroot filesep 'blob_data/sensitivity_test'];

%% 1. parameter setup
% document parameters of interest: 
% -- SST, cloud fraction and surface wind satellite products:
SST_productName = 'gpblend';    % is closer to g16 daily averaged SST than OSTIA.
L4SST_matFN = [SST_productName '_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N'];

% cloud area fraction data (??)
cldfrac_dataFN = [dataroot filesep 'GOES16_L3C_cloudfrac/cloudfrac_5kmres_calculated_from_2km_cloudmask_6-20N_62W-45W.mat'];
% refCF_matFN = 'reference_cloudiness.mat';
% refstate = 'JanFebMean';      %or : 3d_movmean, 5d_movmean


% wind: (both ERA5 and CCMPv2)
ERA5_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/era5_data';
ERA5dataFN = [ERA5_datadir filesep 'era5_WindVectors_DJF2019-2022.nc'];    % I can include full level wind in this file. 

CCMPv2_datadir = './CCMPv2';
CCMPv2_dataFN  = [CCMPv2_datadir filesep 'CCMP_RT_WindVector_Analysis_V2_L3_RSS_JanFeb2020_ATOMIC_broder_dailyMean.nc'];  

% -- specific parameters that can be varied:
% ======== a. feature identification parameters:  ============ %
dT_thres=0.1;                        % temperature anomaly threshould.
dT_thres_arry=[0.1:0.05:0.4];        % use different threshold to get all the blobs out.  we don't need this.
minArea=25;                          % 5x5 pixels: 25km (1pixel is 5km across)
cldcv_thres = 0.1;                   % irrelevant here? threshold to constrain whether or not a warm blob identified is under strong influence of cloud overhead (a.k.a., spatial interpolation.).
visualize_SSTfeatures = false;       % control making movies of detected SST features.

% ======== b. miscelleneous ===== %
search_RadRatio = 3;                 % cutout area size, 2x the bounding box size. for wind speed
feature_norm_type = 'circle';        % or ellipse;

zwind = 10;                          % surface wind is used to rotate the coord.


% ======== c. spatial high pass filtering parameter ============ %
% how to decide this cutoff length scale?? (need some ideas);
HPcutoff_scale_km.SST = 600;
HPcutoff_scale_km.cldfrac = 200;
HPcutoff_scale_km.wind = 200; 

% ========= d. sample random locations or not ========== %
sample_random_blobs = false;
% cldfrac_highpass_cutoff_spatialscale = 600; %[6, 12, 24, 48];
% highpass_cutoff_scale = cldfrac_highpass_cutoff_spatialscale(1);  % 3,  6, 24  % twelve days do not work well with to get the signal out..

% ================ construct caseID from info. above =================== %
xgrd_vec = [-search_RadRatio:0.05:search_RadRatio];              % normalized grid vector range.
ygrd_vec = [-search_RadRatio:0.05:search_RadRatio];

processing_parms = strjoin({['dTthres', num2str(dT_thres)]; ...
    ['minArea', num2str(minArea)]; ...
    %['zwind', num2str(zwind)]; ...
    ['searchRad',num2str(search_RadRatio)]; ...
    %['refCF', refstate]; ...                                   % do i need refCF???
    ['NormType', feature_norm_type]; ...
    ['SSTHPthres', num2str(HPcutoff_scale_km.SST)]; ...
    ['CFHPthres', num2str(HPcutoff_scale_km.cldfrac)]; ...
    ['WindHPthres', num2str(HPcutoff_scale_km.wind)]}, ...
    '_');

% ---- set up caseID and path to store figures and processed data ---- %

caseID = processing_parms;
figsvdir = [casesvdir filesep caseID filesep 'figs'];
datasvdir = [casesvdir filesep caseID filesep 'processed_data'];

if ~exist(figsvdir, 'dir')
    mkdir(figsvdir)
end

if ~exist(datasvdir,'dir')
    mkdir(datasvdir)
end

% ====================   load in all the data:   ======================== %
% --   (ERA5 and CCMPv2 will be loaded inside some functions.
load([dataroot filesep SST_productName filesep L4SST_matFN]);
L4SST_data = data;  clear data

% double check if SST feature movie exist already, if so, change
% visualize_feature flag to false. 
existing_gifs = dir([figsvdir filesep '*.gif']);
egifNames = {existing_gifs.name};

if contains(egifNames, SST_productName)
    % animation file already exist.
    visualize_SSTfeatures = false;
else
    visualize_SSTfeatures = true;
end

%% 2. high-pass filtering (do this first???)


%% 3. feature detection and feature based information extraction, coord. transformation
clear indiv_blobs_collection cloudiness_downwind_FCN SSTcutout_downwind_FCN winddiv_downwind_FCN
if sample_random_blobs
    clear random_blobs cloudiness_downwind_FCN_rb SSTcutout_downwind_FCN_rb winddiv_downwind_FCN_rb
end
for it =1 :length(L4SST_data.time)
    %% a: find daily (averaged) SST warm patches:
    t_now = L4SST_data.time_num(it);
    
    SST0 = L4SST_data.analysed_sst(:,:,it)';
    SSTlon = double(L4SST_data.lon(:,it));
    SSTlat = double(L4SST_data.lat(:,it));
    
    L4res = mean([diff(SSTlon)*111E3.*cosd(SSTlat(1:end-1)), diff(SSTlat)*111E3],'all');
    
    [LON, LAT] = meshgrid(SSTlon, SSTlat);    

    checkflag = false;
    
    % the Guassian filter needs to be double checked with Ju or Brandon.
    % (found a reference that provides a detail description of the FFT
    % based filtering I used.)  
    %% To-do: update the following function to change the function name.
    meanSST_LS =estimate_LSG_SST(SSTlon, SSTlat, SST0, 'method','spectrum','CutoffScale', HPcutoff_scale_km.SST, 'checkflag',checkflag);    % the size is different here.
    
    SST_anom_spatial = SST0 - meanSST_LS;
    fake_cloud_flag = isnan(SST_anom_spatial);
    cloud_flag = fake_cloud_flag;
   
    
    % detect both cold and warm anomalies together.  threshold_type = absolute.
    [blob_info, exitflag] = find_SST_blobs(SST_anom_spatial, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
        'LON_grid', LON, 'LAT_grid', LAT,'checkflag', false, 'threshold_type', 'absolute');
    
    
    %% b: get blob statistics
    %     if thres_cnt > 0
    %         debug = true;
    %         [individual_blobs, Nblobs] = remove_doppelgangers(blob_all, debug);
    %         individual_blobs.time = t_now;
    blob_info.time =t_now;
    blob_info.L4_xres = L4res;
    %   I need to make sure that max_SSTa is correct... (√: fixed)
    %   put in optional varaible to control the cutout size.
    indiv_blobs_collection(it) = retrieve_SST_and_relatedInfo_for_each_warm_blob(blob_info, LON, LAT, SST0, meanSST_LS, search_RadRatio);
    
    if sample_random_blobs
        %--- get random circular regions of interest (ROIs) that have the same size as
        %--- the realistic blobs in parallel.
        random_blobs(it) = generate_random_bloblocs_with_realistic_blobstats(SSTlon, SSTlat, indiv_blobs_collection(it));
    end
    % visualize it better and save as animations:
    if visualize_SSTfeatures
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
    
    
    % note: I also need to do the coordinate transformation for L4 SST
    % features.
    %% c: get averaged ERA5 wind direction within 2x bounding box area. (need a utility function here)
    %%    we could compute the surface wind divergence as well in the bounding box region. (need a utility function here)
    % TO-DO : update function to allow specification of the level of
    % ERA5 wind to use.  (done √)
    % use ERA5 for large scale (daily averaged) wind information.
    % is the wind information averaged over the SST cutouts area?
    [ave_windInfo, SearchArea_dw] = get_ERA5_windInfo_over_features(t_now, indiv_blobs_collection(it), search_RadRatio, ERA5dataFN, zwind);      % output: structure: dir, spd, u,v;
    %[ave_windInfo, SearchArea_dw] = get_CCMPv2_windInfo_over_features(t_now, indiv_blobs_collection(it), search_RadRatio, ERA5dataFN, zwind);      % output: structure: dir, spd, u,v;
    
    % using ave_windInfo and the daily SST map, I can compute effective
    % SST gradient:  % the effective SST gradient is computed using the
    % background averaged wind speed...
    effective_SSTgrad = compute_effective_SSTgrad_over_features(ave_windInfo, indiv_blobs_collection(it));
    
    % add compute laplacian using the same method as used in Li and Carbone (2012)
    % laplacianSST = compute_laplacian_SST_over_features(ave_windInfo, indiv_blobs_collection(it))
    
    %% d: get L3C cloud mask and compute cloudiness for the area within 2x bounding box. (need a utility function here)
    % TO-DO: add enhanced cloudiness from this function as well. (done √)
    %cloudiness = compute_cloudiness_over_features_from_L3C_cloudmasks(t_now, SearchArea_dw, indiv_blobs_collection(it), refCF_matFN, refstate);  % cloudiness is a structure, contain cloudiness over each bounding box.
    
    % this cloud frac is only for the 8~18N, -59~-48W...

    cloudfrac = compute_cloudfrac_over_features_from_L3C_cloudmasks(cldfrac_dataFN, t_now, SearchArea_dw, indiv_blobs_collection(it), HPcutoff_scale_km.cldfrac);

    %% e: get wind aligned, feature centric, normalized coordinate to show the cloudiness map. (need a utility function here)
    % TO-DO: add option to approx feature by circle instead of by
    % ellipse. (done √)
    
    % To-DO: generalize the function to map input data (cloudiness or
    % SST) to new Coord.
    cloudiness_downwind_FCN = map_data_to_newCoord(ave_windInfo,  cloudfrac, indiv_blobs_collection(it), 'options', feature_norm_type, ...
        'xgrid', xgrd_vec, 'ygrid',ygrd_vec);
    

    % construct a variable that contain SST cutouts for each blob.
    SSTfieldN = {'LON_cutouts', 'LAT_cutouts', 'SST_cutouts', 'SSTa_cutouts'};
    for ff = 1:length(SSTfieldN)
        FN = SSTfieldN{ff};
        nblobs = length(indiv_blobs_collection(it).(FN));
        for ib = 1:nblobs
            SSTblob_cutouts(ib).(FN) = indiv_blobs_collection(it).(FN){ib};
        end
    end
    for ib = 1:nblobs
        SSTblob_cutouts(ib).eSSTgrad = effective_SSTgrad(ib).val;
    end
    
    % I could increase the cutout size (the same cutout radius as the
    % cloudiness?
    SSTcutout_downwind_FCN = map_data_to_newCoord(ave_windInfo,   SSTblob_cutouts, indiv_blobs_collection(it), 'options', feature_norm_type, ...
        'xgrid', xgrd_vec, 'ygrid', ygrd_vec);
    
    
    %% f. compute the ERA5 wind convergence in the normalized coordinate. (can be done later)
    winddiv = compute_winddiv_over_features(CCMPv2_dataFN, t_now, SearchArea_dw, indiv_blobs_collection(it), HPcutoff_scale_km.wind);
    winddiv_downwind_FCN = map_data_to_newCoord(ave_windInfo,  winddiv, indiv_blobs_collection(it), 'options', feature_norm_type, ...
        'xgrid', xgrd_vec, 'ygrid',ygrd_vec);
    
    clear cloudfrac SSTblob_cutouts winddiv
    
    
    if sample_random_blobs
        [ave_windInfo_rb, SearchArea_dw_rb] = get_ERA5_windInfo_over_features(t_now, random_blobs(it), search_RadRatio, ERA5dataFN, zwind);      % output: structure: dir, spd, u,v;
        
        effective_SSTgrad_rb = compute_effective_SSTgrad_over_features(ave_windInfo_rb, random_blobs(it));
        
        cloudfrac_rb = compute_cloudfrac_over_features_from_L3C_cloudmasks(cldfrac_dataFN, t_now, SearchArea_dw_rb, random_blobs(it), HPcutoff_scale_km.cldfrac);
        
        cloudiness_downwind_FCN_rb = map_data_to_newCoord(ave_windInfo_rb,  cloudfrac_rb, random_blobs(it), 'options', feature_norm_type, ...
            'xgrid', xgrd_vec, 'ygrid',ygrd_vec);
        
        
        SSTfieldN = {'LON_cutouts', 'LAT_cutouts', 'SST_cutouts', 'SSTa_cutouts'};
        for ff = 1:length(SSTfieldN)
            FN = SSTfieldN{ff};
            nblobs = length(random_blobs(it).(FN));
            for ib = 1:nblobs
                SSTblob_cutouts_rb(ib).(FN) = random_blobs(it).(FN){ib};
            end
        end
        for ib = 1:nblobs
            SSTblob_cutouts_rb(ib).eSSTgrad = effective_SSTgrad_rb(ib).val;
        end
        
        
        SSTcutout_downwind_FCN_rb = map_data_to_newCoord(ave_windInfo_rb,   SSTblob_cutouts_rb, random_blobs(it), 'options', feature_norm_type, ...
            'xgrid', xgrd_vec, 'ygrid', ygrd_vec);
        
        winddiv_rb = compute_winddiv_over_features(CCMPv2_dataFN, t_now, SearchArea_dw_rb, random_blobs(it), HPcutoff_scale_km.wind);
        winddiv_downwind_FCN_rb = map_data_to_newCoord(ave_windInfo_rb,  winddiv_rb, random_blobs(it), 'options', feature_norm_type, ...
            'xgrid', xgrd_vec, 'ygrid',ygrd_vec);
        
        clear cloudfrac_rb SSTblob_cutouts_rb winddiv_rb
    end
    
    
    close all
    
    %% the following data are needed for composite purposes. 
    % rearrange data for storaging purpose
    downwind_FCN_regridded(it).XX = SSTcutout_downwind_FCN(1).regridded.XX;
    downwind_FCN_regridded(it).YY = SSTcutout_downwind_FCN(1).regridded.YY; 
    
    downwind_FCN_regridded(it).ave_windInfo = reorganize_struct_to_array(ave_windInfo);
    downwind_FCN_regridded(it).SSTInfo = reorganize_struct_to_array([SSTcutout_downwind_FCN.regridded]);   % can I change this structure into a matlab array instead?
    downwind_FCN_regridded(it).cldfrac = reorganize_struct_to_array([cloudiness_downwind_FCN.regridded]);
    downwind_FCN_regridded(it).winddiv = reorganize_struct_to_array([winddiv_downwind_FCN.regridded]);
    
    if sample_random_blobs
        % for the random blob.
        rb_downwind_FCN_regridded(it).XX = SSTcutout_downwind_FCN_rb(1).regridded.XX;
        rb_downwind_FCN_regridded(it).YY = SSTcutout_downwind_FCN_rb(1).regridded.YY;
        rb_downwind_FCN_regridded(it).ave_windInfo = reorganize_struct_to_array(ave_windInfo_rb);
        rb_downwind_FCN_regridded(it).SSTInfo = reorganize_struct_to_array([SSTcutout_downwind_FCN_rb.regridded]);
        rb_downwind_FCN_regridded(it).cldfrac = reorganize_struct_to_array([cloudiness_downwind_FCN_rb.regridded]);
        rb_downwind_FCN_regridded(it).winddiv = reorganize_struct_to_array([winddiv_downwind_FCN_rb.regridded]);
    end
    %% other things that I also need: information at the physical space;  (plot case by case)
    % remove the cell array. 
    tmp = rmfield(SSTcutout_downwind_FCN,{'WindAligned_Coord', 'mean_wspd','regridded'});
    cutouts_physicalspace(it).SSTInfo = tmp;
    
    tmp = rmfield(cloudiness_downwind_FCN,{'WindAligned_Coord', 'mean_wspd','regridded'});   
    cutouts_physicalspace(it).cldfrac = tmp;
    
    tmp = rmfield(winddiv_downwind_FCN,{'WindAligned_Coord', 'mean_wspd','regridded'});       
    cutouts_physicalspace(it).winddiv = tmp;
    

end


%% 4. feature based composite 
run('feature_based_composite_analysis_plugin_version.m');

%% 5. data and result storage
% try to see what I can combine together to save storage space.
% --------------------------------------------------------------------- %
% ==== save animations ==== %
if visualize_SSTfeatures
    if exist('images', 'var')
        gifname = [SST_productName 'SST_warm_and_cold_features.gif'];
        gif_abspath = [figsvdir filesep gifname];
        write_frames_into_video(images, gif_abspath, 'gif');
    end
end

% saving in a different format.
% ==== save data ==== % 
% everything in just 1 matlab data file. 
newformat = true;
OutPutFN = [SST_productName '_SST-features_properties.mat'];
save([datasvdir filesep OutPutFN], 'indiv_blobs_collection', 'downwind_FCN_regridded', 'cutouts_physicalspace', '-v7.3');

if sample_random_blobs
    OutPutFN_rb = [SST_productName 'random_ROIs_for_siglev_estimation.mat'];
    save([datasvdir filesep OutPutFN_rb], 'random_blobs', 'rb_downwind_FCN_regridded','-v7.3');
end

%% the following is the old format 
if ~newformat
    dataFN_prefix = [SST_productName '_daily_blobstats'];
    dataFN = [strjoin({dataFN_prefix, processing_parms}, '_') '.mat'];
    save([datasvdir filesep dataFN],'indiv_blobs_collection','-v7.3');
    
    dataFN_prefix = [SST_productName '_daily_SSTcutouts'];
    dataFN = [strjoin({dataFN_prefix, processing_parms}, '_') '.mat'];
    save([datasvdir filesep dataFN],'SSTcutout_downwind_FCN','-v7.3');
    
    dataFN_prefix = 'GOES16_L3C_5kmCldFrac';
    dataFN = [strjoin({dataFN_prefix, processing_parms}, '_') '.mat'];
    save([datasvdir filesep dataFN],'cloudiness_downwind_FCN','-v7.3');
    
    dataFN_prefix = 'CCMPv2Wind';
    dataFN = [strjoin({dataFN_prefix, processing_parms}, '_') '.mat'];
    save([datasvdir filesep dataFN],'winddiv_downwind_FCN','-v7.3');
end

% == what else? == %