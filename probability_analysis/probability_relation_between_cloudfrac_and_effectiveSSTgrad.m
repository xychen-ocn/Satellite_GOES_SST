% Description:
% Authors:
% Dates:


clear all; clc; 

%% set up environment:
% add paths:
addpath('/Users/xchen/Documents/GitHub/SAM_LES/matlab');
dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
addpath(genpath([dataroot filesep 'util']));              
addpath([dataroot filesep 'probability_analysis']);

%% I/O: 
% define file names to be used:
dataFolder.GOES = 'gpblend'; 
dataFolder.ERA5 = 'era5_data';

ERA5_dataFN = 'era5_10mWindVectors_DJF2019-2022.nc';
L4_dataFN = 'gpblend_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N.mat';

% -- L4 SST data:
load([dataroot filesep dataFolder.GOES filesep L4_dataFN]);
L4data = data;  clear data

% -- ERA5 data:
absFN = [dataroot filesep dataFolder.ERA5 filesep ERA5_dataFN];
ERA5data=read_netCDF_into_matlab_structure(absFN);


% select subregion:
subregion = [-59, -48; 8 18];

%% data processing:
% -- L4 SST data:
%    compute SST gradient from L4 data: 
cutoff_scale_km = 600;
for it = 1:length(L4data.time)
    SST0 = L4data.analysed_sst(:,:,it)';
    lon = double(L4data.lon(:,it));
    lat = double(L4data.lat(:,it));
    
    midlat = (lat(1:end-1) + lat(2:end))/2;
    
    xvec = [0; cumsum(diff(lon)*111E3.*cosd(midlat))];   % km    % I am doubting this ... double check.
    yvec = [0; cumsum(diff(lat)*111E3)];                 % km
    
    meanSST_LS(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',false);    % the size is different here.
    
    SST_anom_spatial = SST0 - meanSST_LS;
%     SSTgrad(it) = spatial_grad(SST_anom_spatial, xvec, yvec);  %units: K/km
%     SSTLaplacian(:,:, it) = spatial_laplacian(SST_anom_spatial, xvec, yvec);
    
    SSTgrad(it) = spatial_grad(SST0, xvec, yvec);  %units: K/km
    SSTLaplacian(:,:, it) = spatial_laplacian(SST0, xvec, yvec);
    
end

% build two meshes at L4 resolution.
[LON_L4, LAT_L4] = meshgrid(L4data.lon(:,1), L4data.lat(:,1));
lonmask = L4data.lon(:,1)>=subregion(1,1) & L4data.lon(:,1)<=subregion(1,2);
latmask = L4data.lat(:,1)>=subregion(2,1) & L4data.lat(:,1)<=subregion(2,2);

LON_L4sub = LON_L4(latmask, lonmask); 
LAT_L4sub = LAT_L4(latmask, lonmask); 

% -- used ERA5 wind to compute effective SST gradient.
ERA5data.wdir = atan2(ERA5data.vwnd, ERA5data.uwnd);

% compute uvec * SST grad 
% deal with resolution: from coarse ERA5 grid to fine L4 grid. 
[ERA5data.LON, ERA5data.LAT] = meshgrid(ERA5data.lon, ERA5data.lat);
for it = 1:length(ERA5data.time)
    ERA5_uc_L4(:,:,it) = interp2(ERA5data.LON-360, ERA5data.LAT, ERA5data.uwnd(:,:,it)', LON_L4, LAT_L4);
    ERA5_vc_L4(:,:,it) = interp2(ERA5data.LON-360, ERA5data.LAT, ERA5data.vwnd(:,:,it)', LON_L4, LAT_L4);
end

dx = diff(LON_L4,1, 2);
dx = [zeros(size(dx,1),1)  dx];

dy = diff(LAT_L4,1,1);
dy = [zeros(1,size(dy,2)); dy];

XX_L4 = 111*cosd(LAT_L4).*cumsum(dx,2);     % km
YY_L4 = 111.*cumsum(dy,1);                        % km

for it = 1:length(ERA5data.time)
    ERA5_wdiv_L4(:,:,it) = divergence(XX_L4, YY_L4, ERA5_uc_L4(:,:,it), ERA5_vc_L4(:,:,it));
end


for it = 1:length(L4data.time_num)
    tid_era5 = find(ERA5data.time_num==floor(L4data.time_num(it)));
    air_Ttend(:,:,it) = (ERA5_uc_L4(:,:,tid_era5).*SSTgrad(it).xcomp + ERA5_vc_L4(:,:,tid_era5) .* SSTgrad(it).ycomp)*3.6;  % units: K/h; (1m/s = 3.6km/h)
    winddiv(:,:,it) = ERA5_wdiv_L4(:,:,tid_era5);
end


% -- cloud mask data: I can get cloudiness and cloud fraction from the
% cloud mask data
L3C_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/g16';
% load in L3C daily dataset;
for it = 1:length(L4data.time)
    t = L4data.time_num(it);
    tday = floor(t);
    dateID = datestr(tday, 'mmmdd');
    DOY = tday - datenum(2019, 12, 31);
    
    % daily cloud mask.
    matFN = ['GOES_SST_' num2str(DOY,'%3.3i') '_' dateID '_lon62W_TO_42W_lat5N_TO_25N.mat'];
    disp(['loading' dateID])
    load([L3C_datadir filesep matFN]);
    
    [GOES_LON, GOES_LAT] = meshgrid(GOES_ATOMIC.lon(:,1), GOES_ATOMIC.lat(:,1));
    lonmask = GOES_ATOMIC.lon(:,1)>=subregion(1,1) & GOES_ATOMIC.lon(:,1)<=subregion(1,2);
    latmask = GOES_ATOMIC.lat(:,1)>=subregion(2,1) & GOES_ATOMIC.lat(:,1)<=subregion(2,2);

    SST_all = permute(GOES_ATOMIC.sea_surface_temperature, [2,1,3]);   % dimension is Nlat, Nlon, ..
    time_all = GOES_ATOMIC.time_num;
    
    % check cloud mask and throw away maps that contains all cloudy pixels.
    %% QC data:
    % throw away maps that is 100% garbage;
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
    
    if length(time_all)>=12        % half of the day has valid cloud mask.
        cldfreq(:,:,it) = compute_cloudfreq(SST_all);
    else
        cldfreq(:,:,it) = nan(size(SST_all));
    end
    
    
    %% Note: this step needs to be further checked; Instead of cloud frequency , perhaps cloud coverage should be used to describe cloudiness...
    % select the area of interest to reduce grid point .
    native_grid.LON = GOES_LON(latmask, lonmask);
    native_grid.LAT = GOES_LAT(latmask, lonmask);
    coarse_grid.LON = LON_L4sub;
    coarse_grid.LAT = LAT_L4sub;
    native_cldmask = double(isnan(SST_all(latmask, lonmask)));
    
    cldfrac_tmp = compute_spatial_cloud_fraction_in_coarse_grid(native_grid, native_cldmask, coarse_grid);
    cloudiness_ds(it).cldfrac_hrly = cldfrac_tmp;
    cloudiness_ds(it).cldfrac_dailymean = mean(cldfrac_tmp, 3,'omitnan');
end

% n0 = 1;
% for it = 1:length(cldfrac)
%     nn = n0-1 + size(cldfrac{it},3);
%     cldfrac_2mon(:,:,n0:nn) = cldfrac{it};
%     n0 = nn+1;
% end
cldfrac_2mon = [cloudiness_ds.cldfrac_hrly];

% find averaged cloud fraction: (spatial + temporal)
ave_cloudfrac.spatemp_2months = mean(cldfrac_2mon,'all', 'omitnan');
ave_cloudfrac.temp_2months = mean(cldfrac_2mon, 3, 'omitnan');             % a function of space


% now look at the distribution according to different large scale
% conditions:
dayType = get_dayType_by_U10_and_LTS;
cloudtypes = fieldnames(dayType);
for tt = 1:4
    CN = cloudtypes{tt};
    tids = ismember(floor(L4data.time_num), dayType.(CN));
    aTtend.(CN) = air_Ttend(latmask, lonmask, tids);
    %CFc.(CN) = cldfreq_coarse(latmask, lonmask, tids);
    tmp = cloudiness_ds(tids);
    CldFrac.(CN) = [tmp.cldfrac_dailymean];
    SSTlap.(CN) = SSTLaplacian(latmask, lonmask, tids);
    wdiv.(CN) = winddiv(latmask, lonmask, tids);
    SST_LS.(CN) = meanSST_LS(latmask,lonmask, tids);
    
end

%% make plots:
% relative change or anomalies in cloud fraction in different atmospheric
% conditions:
% ----------------------------------------------------------------------- %
binwidth = 5;
xPrcTileBins = 0:binwidth:100;
xbincen = 0.5*(xPrcTileBins(1:end-1) + xPrcTileBins(2:end));

% -- 1. for all data:
cldfrac_all = [cloudiness_ds.cldfrac_dailymean];
[aRCC_l, aRCC_dm, xedge] = compute_relative_change_of_cloudiness_in_xbins(air_Ttend(latmask, lonmask,:), xPrcTileBins, cldfrac_all, ave_cloudfrac.spatemp_2months, LON_L4sub, LAT_L4sub); 

% do an estimation of 5% null hypothesis level:
% randomly took 5% of data 1000 times, and find the standard deviation of
% average cloudiness; use 1.96*stdv as the range.
null_all= bootstrap_nullhypothesis_level_for_relative_cloudiness_change(air_Ttend(latmask, lonmask,:), binwidth, cldfrac_all, ave_cloudfrac.spatemp_2months);

    
prcs_1 = prctile(air_Ttend(latmask, lonmask,:), xPrcTileBins,'all');
prc_rank = interp1(prcs_1, xPrcTileBins, 0);

figure(11); clf;
hold on
bar(xbincen, aRCC_l,1.0, 'FaceAlpha', 0.5);
plot([prc_rank, prc_rank],[-15, 30], '--b', 'linewidth', 1.2);
hold on;
% plot 95% confidence level of the mean change of CF in 5% of data; 
vlocs = [0, null_all.stdv*1.96 + null_all.mean; 
         100, null_all.stdv*1.96 + null_all.mean; 
         100, -null_all.stdv*1.96 + null_all.mean; 
         0, -null_all.stdv*1.96 + null_all.mean; ];
         
patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');    
xlim([0 100]);    
xlabel('\bf{u} \cdot \nablaSST percentiles');
ylabel('average change of cloudiness in each percentiles (%)');
set(gca,'fontsize',14);
title('Jan and Feb 2020 (all)');
grid on

ax1 = gca;
ax2 = axes('pos', ax1.Position);
bar(ax2, xbin, RCC_l,1.0, 'FaceAlpha', 0);

ax2.XAxisLocation = 'top';
ax2.YColor='none';
ax2.Color= 'none';
ax2.XLim = ax1.XLim;

for i = 1:length(ax2.XTick)
    xticklabels{i} = num2str(prctile(air_Ttend(latmask, lonmask,:), ax2.XTick(i),'all'), '%.2e');
end
ax2.XTickLabel =xticklabels;
set(ax2, 'fontsize',14)
    
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_CloudFrac_vs_effective_SSTgrad_percentiles.jpg', [0 0 10 8]); 


% -- 2. for each atmospheric regime:

for tt = 1:4
    CN = cloudtypes{tt};
    [RCC_l(tt), RCC_dm(tt)] = compute_relative_change_of_cloudiness_in_xbins(aTtend.(CN), xPrcTileBins, CFc.(CN), ave_cloudfrac.spatemp_2months, LON_L4sub, LAT_L4sub);
    null(tt) = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(aTtend.(CN), binwidth, CFc.(CN), ave_cloudfrac.spatemp_2months, 1000);
end

figure(10); clf
for tt = 1:4
    CN = cloudtypes{tt};
    prcs_1 = prctile(aTtend.(CN), xPrcTileBins,'all');
    prc_rank = interp1(prcs_1, xPrcTileBins, 0);
    
    subplot(2,2,tt);
    hold on
    %bar(xbin, RCC_dm, 'FaceAlpha',1);
    bar(xbincen, RCC_l(tt), 1.0,'FaceAlpha', 0.5);
    % plot 95% confidence level of the mean change of CF in 5% of data;
    vlocs = [0, null(tt).stdv*1.96 + null(tt).mean;
        100, null(tt).stdv*1.96 + null(tt).mean;
        100, -null(tt).stdv*1.96 + null(tt).mean;
        0, -null(tt).stdv*1.96 + null(tt).mean; ];
    
    patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');


    xlim([0 100]);
    xlabel('\bf{u} \cdot \nablaSST percentiles');
    ylabel('average change of cloudiness (%)');
    set(gca,'fontsize',14);
    %title(CN);
    grid on
    yrange = get(gca,'ylim');
    plot([prc_rank, prc_rank],yrange, '--k');
    hold on;
    
    ax1 = gca;
    ax2 = axes('pos', ax1.Position);
    bar(ax2, xbin, RCC_l,1.0, 'FaceAlpha', 0);

    ax2.XAxisLocation = 'top';
    ax2.YColor='none';
    ax2.Color= 'none';
    ax2.XLim = ax1.XLim;
    ax2.YLim = ax1.YLim;
    
    for i = 1:length(ax2.XTick)
        xticklabels{i} = num2str(prctile(aTtend.(CN), ax2.XTick(i),'all'), '%.2e');
    end
    ax2.XTickLabel =xticklabels;
    ax2.XTickLabelRotation=25;
    set(ax2, 'fontsize',14)

end
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_CloudFrac_vs_effective_SSTgrad_percentiles_4types_v2.jpg', [0 0 10 8]); 
