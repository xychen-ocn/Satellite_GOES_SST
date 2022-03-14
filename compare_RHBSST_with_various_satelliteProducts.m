% This script is used to compare various satellite SST product with the RHB
% sampled data. The results will be displayed using Taylor diagram.
% Goal: provide statistical measures to show that
% GOES-16 SST captures the RHB sampled variability of SST well.
%
% input RHB data: segments on days where the RHB sampled high SST variance
% (>75th prcentile)
%
% date: Jan 26, 2022; (drafted); need advise from Gary Wick to do the
% analysis in the right direction to move further.

%% load in the RHB days of interests for comparisons:
RHB_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
load([RHB_datadir filesep 'EUREC4A_ATOMIC_RHB_AllinOne_v1.3.mat']);

strong_SSTvar_dates = rhbdates.strong_SSTvar;
moderate_SSTvar_dates = rhbdates.moderate_SSTvar;

basetimenum = datenum(2019,12,31);
all_DOYs = [strong_SSTvar_dates, moderate_SSTvar_dates]-basetimenum;
all_DOYs = sort(all_DOYs);



%% collocate SST measurements with RHB sampled SST
for it = 1:length(all_DOYs)
    % load in satellite data:
    DOY = all_DOYs(it);
    dateID = datestr(basetimenum + (DOY),'mmmdd');
    load(['GOES_SST_' dateID '.mat']);
    disp(['working on ' dateID]);
    
    QC_mask = ones(size(GOES_ATOMIC.quality_level));
    QC_mask(GOES_ATOMIC.quality_level~=5) = NaN;
    masked_SST = GOES_ATOMIC.sea_surface_temperature.*QC_mask;
    
    time_mask =(RHB.time>=basetimenum+(DOY)) & (RHB.time<basetimenum+DOY+1);
    RHB_lon = RHB.lon(time_mask);
    RHB_lat = RHB.lat(time_mask);
    
    RHB_time = RHB.time(time_mask);
    
    
    % so now, we can do some interpolation and compare the SST:
    lon_sat = double(GOES_ATOMIC.lon(:,1));
    lat_sat = double(GOES_ATOMIC.lat(:,1));
    time_sat = GOES_ATOMIC.time_num;
    SST_sat = permute(GOES_ATOMIC.sea_surface_temperature,[2,1,3]);
    
    
    satSST_atRHB = interp3(lon_sat, lat_sat,time_sat, SST_sat, RHB_lon, RHB_lat, RHB_time);
    SST_collocated.(dateID).G16 = satSST_atRHB;
    SST_collocated.(dateID).RHB =  RHB.tsea(time_mask);
    
    % compute statistics:
    stats.G16_stdv(it) = std(satSST_atRHB,1, 'omitnan');
    stats.RHB_stdv(it) = std(RHB.tsea(time_mask),1,'omitnan');
    stats.CorrCoef(it) = get_CorrCoef(satSST_atRHB, RHB.tsea(time_mask));
    stats.cRMSE(it) = get_centered_RMSE(satSST_atRHB, RHB.tsea(time_mask));
    stats.RMSE(it) = get_RMSE(satSST_atRHB, RHB.tsea(time_mask));
end

figure
subplot(1,2,1)
plot(stats.G16_stdv, stats.RHB_stdv,'.k','markerSize',12);
hold on
plot([0.:0.1:0.3],[0.:0.1:0.3],'-b');
hold on
xlabel('GOES-16');
ylabel('RHB');

subplot(1,2,2)
plot(stats.CorrCoef, stats.cRMSE,'.k','markerSize',12)
xlabel('Correlation');
ylabel('centered RMSE');


%% compute standard deviation for both the RHB and the satellite, compute correlation, and centered RMSE.
stds = [stats.RHB_stdv, stats.G16_stdv];
%RMSE = [0, stats.cRMSE];
Corr = [1, stats.CorrCoef];

RMSE_should_be = sqrt(stds.^2 + stds(1)^2 -2*stds.*stds(1).*Corr);

% what is the better way to show that the Satellite data can capture well
% the SST variability in this region? (ask Gary).


% stds = [0.3, 0.5];
% RMSE = [0, 0.3391];
% Corr = [1. 0.75];
taylordiag(stds, RMSE_should_be, Corr);
