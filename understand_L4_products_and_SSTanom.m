% Purpose: find warm anomalies from L4 data;
clear all; clc; close all;

dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';

matFN = 'TemporalAveraged_LowPassFiltered_thres600km_SST.mat';
matFN2 = '_L4_daily_JanFeb_ATOMIC_lon62W_TO_45W_lat5N_TO_20N';

gpbdata = [dataroot filesep 'gpblend' filesep matFN];
ostdata = [dataroot filesep 'ostia' filesep matFN];
g16data = [dataroot filesep 'g16' filesep matFN];

% load daily maps to compare the daily SST maps:
gpbdata2 = [dataroot filesep 'gpblend' filesep 'gpblend' matFN2];
ostdata2 = [dataroot filesep 'ostia' filesep 'ostia' matFN2];
g16data2 = [dataroot filesep 'g16' filesep 'smaller_domain/G16_daily_averaged_SST.mat'];

% compare how different is the mean field (leftover question from previous sensitivity investigation)
load(gpbdata);
load(g16data);
load(ostdata);

[LONG16, LATG16] = meshgrid(G16_movaved.lon, G16_movaved.lat);
[LONL4, LATL4] = meshgrid(GPBLEND_movaved.lon, GPBLEND_movaved.lat);

for it = 1:size(G16_movaved.wndsz_3d.SST,3)    
    G16_coarsed.wndsz_3d(:,:,it) = interp2(LONG16, LATG16, squeeze(G16_movaved.wndsz_3d.SST(:,:,it))', LONL4, LATL4);
end

for it = 1:size(G16_movaved.wndsz_5d.SST,3)    
    G16_coarsed.wndsz_5d(:,:,it) = interp2(LONG16, LATG16, squeeze(G16_movaved.wndsz_5d.SST(:,:,it))', LONL4, LATL4);
end

%% daily maps:
gpb = load(gpbdata2);
ost = load(ostdata2);
load(g16data2);          % daily averages:

% fine G16 to coarse grids for comparison:
[LONG16s, LATG16s] = meshgrid(G16_daily_ave.lon(:,1), G16_daily_ave.lat(:,1));
[LONL4s, LATL4s] = meshgrid(gpb.data.lon(:,1), gpb.data.lat(:,1));

for it = 1: length(G16_daily_ave.time)
    G16_coarsed.daily(:,:,it) = interp2(LONG16s, LATG16s, squeeze(G16_daily_ave.SST(:,:,it))', LONL4s, LATL4s);
end


% compute the averaged mean difference over the two month;
% compute the variance of the mean difference:
wz = 1;
if wz>1
    WN = ['wndsz_' num2str(wz) 'd'];
    
    
    for it = 1:size(G16_movaved.(WN).SST,3)
        G16_coarsed.(WN)(:,:,it) = interp2(LONG16, LATG16, squeeze(G16_movaved.(WN).SST(:,:,it))', LONL4, LATL4);
    end
    
    tmask = ismember(G16_movaved.(WN).time, floor(GPBLEND_movaved.(WN).time));
    G16_GPB_diff = permute(GPBLEND_movaved.(WN).SST,[2,1,3])  - G16_coarsed.(WN)(:,:,tmask);
    
    tmask = ismember(G16_movaved.(WN).time, floor(OSTIA_movaved.(WN).time));
    G16_OSTIA_diff = permute(OSTIA_movaved.(WN).SST, [2,1,3]) - G16_coarsed.(WN)(:,:,tmask);
    
    tmask = ismember(OSTIA_movaved.(WN).time, GPBLEND_movaved.(WN).time);
    OSTIA_GPB_diff = permute(GPBLEND_movaved.(WN).SST - OSTIA_movaved.(WN).SST(:,:,tmask),[2,1,3]);
    
    LON = LONL4; LAT = LATL4;
elseif wz==1
    WN = 'daily';
    
    tmask = ismember(G16_daily_ave.time, floor(gpb.data.time_num));
    G16_GPB_diff = permute(gpb.data.analysed_sst, [2, 1,3]) - G16_coarsed.(WN)(:,:,tmask);
    
    tmask = ismember(G16_daily_ave.time, floor(ost.data.time_num));
    G16_OSTIA_diff = permute(ost.data.analysed_sst, [2,1,3]) - G16_coarsed.(WN)(:,:,tmask);
    
    tmask = ismember(ost.data.time_num, gpb.data.time_num);
    OSTIA_GPB_diff = permute(gpb.data.analysed_sst - ost.data.analysed_sst(:,:,tmask),[2,1,3]);
    
    LON = LONL4s; LAT = LATL4s;

end

mean_diff.GPBfromG16 = mean(G16_GPB_diff, 3, 'omitnan');
mean_diff.OSTIAfromG16 = mean(G16_OSTIA_diff, 3, 'omitnan');
mean_diff.OSTIAfromGPB = mean(OSTIA_GPB_diff,3, 'omitnan');

% compute the spatial variance of the difference:
var_diff.GPBfromG16 = squeeze(var(G16_GPB_diff,1,3, 'omitnan')); 
var_diff.OSTIAfromG16 = squeeze(var(G16_OSTIA_diff, 1, 3, 'omitnan'));
var_diff.OSTIAfromGPB = squeeze(var(OSTIA_GPB_diff, 1, 3, 'omitnan'));

var_diff_ts.GPBfromG16 = squeeze(var(G16_GPB_diff,1,[1,2], 'omitnan')); 
var_diff_ts.OSTIAfromG16 = squeeze(var(G16_OSTIA_diff, 1, [1,2], 'omitnan'));
var_diff_ts.OSTIAfromGPB = squeeze(var(OSTIA_GPB_diff, 1, [1,2], 'omitnan'));


mean_diff_ts.GPBfromG16 = squeeze(mean(G16_GPB_diff,[1,2], 'omitnan'));
mean_diff_ts.OSTIAfromG16 = squeeze(mean(G16_OSTIA_diff,[1,2], 'omitnan'));
mean_diff_ts.OSTIAfromGPB = squeeze(mean(OSTIA_GPB_diff,[1,2], 'omitnan'));


% compare 3d and 5d averaged mean field:
%% checking differences for 1 day: small differences between GPBLEND and OSTIA;
%  GPBlend and OSTIA both slightly cooler than the G16 (~0.5K), 
%  variance is smaller for GPB overall. and cool less in the lower left
%  corner of the domain.
if 1==0
    
    figure(1);
    subplot(3,3,1)
    pcolor(G16_movaved.lon, G16_movaved.lat, G16_movaved.wndsz_3d.SST(:,:,1)'-273.15);
    shading flat;
    axis([-59, -48, 8 18]);
    colorbar
    caxis([25, 29]);
    title('G16');
    
    subplot(3,3,2)
    pcolor(GPBLEND_movaved.lon, GPBLEND_movaved.lat, GPBLEND_movaved.wndsz_3d.SST(:,:,1)'-273.15);
    shading flat;
    axis([-59, -48, 8 18]);
    colorbar
    caxis([25, 29]);
    title('GPBlend');
    
    
    hsub = subplot(3,3,3);
    pcolor(GPBLEND_movaved.lon, GPBLEND_movaved.lat, GPBLEND_movaved.wndsz_3d.SST(:,:,1)'- G16_coarsed_wndsz_3d(:,:,1));
    shading flat;
    axis([-59, -48, 8 18]);
    caxis([-1 1]);
    colormap(hsub, redblue);
    colorbar(hsub);
    title('GPBlend - G16');
    
    
    
    subplot(3,3,5)
    pcolor(OSTIA_movaved.lon, OSTIA_movaved.lat, OSTIA_movaved.wndsz_3d.SST(:,:,1)'-273.15);
    shading flat;
    axis([-59, -48, 8 18]);
    colorbar
    caxis([25, 29]);
    title('OSTIA');
    
    
    hsub =subplot(3,3,6);
    pcolor(GPBLEND_movaved.lon, GPBLEND_movaved.lat, OSTIA_movaved.wndsz_3d.SST(:,:,1)'- G16_coarsed_wndsz_3d(:,:,1));
    shading flat;
    axis([-59, -48, 8 18]);
    caxis([-1 1]);
    colormap(hsub, redblue);
    colorbar(hsub);
    title('OSTIA - G16');
    
    hsub = subplot(3,3,8);
    pcolor(GPBLEND_movaved.lon, GPBLEND_movaved.lat, GPBLEND_movaved.wndsz_3d.SST(:,:,1)'- OSTIA_movaved.wndsz_3d.SST(:,:,1)');
    shading flat;
    axis([-59, -48, 8 18]);
    caxis([-1 1]);
    title('GPBLEND - OSTIA');
    colormap(hsub, redblue);
    colorbar(hsub);
    
end

%% check the averaged differences in these two months;
figure(2); clf; 
hsub = subplot(3,3,1);
pcolor(LON, LAT, mean_diff.GPBfromG16); shading flat;
shading flat; 
caxis([-0.5 0.5]);
title('GPBLEND - G16');
colormap(hsub, redblue);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);

hsub = subplot(3,3,2);
pcolor(LON, LAT, mean_diff.OSTIAfromG16); shading flat;
shading flat; 
caxis([-0.5 0.5]);
title('OSTIA - G16');
colormap(hsub, redblue);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);


hsub = subplot(3,3,3);
pcolor(LON, LAT, mean_diff.OSTIAfromGPB); shading flat;
shading flat; 
axis([-59, -48, 8 18]);
caxis([-0.5 0.5]);
title('GPBLEND - OSTIA');
colormap(hsub, redblue);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);

% plot local temporal variance:
hsub = subplot(3,3,4);
pcolor(LON, LAT, var_diff.GPBfromG16); shading flat;
shading flat; 
caxis([0 0.055]);
title('variance of GPBLEND-G16');
colormap(hsub, parula);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);

hsub = subplot(3,3,5);
pcolor(LON, LAT, var_diff.OSTIAfromG16); shading flat;
shading flat; 
caxis([0 0.055]);
title('variane: OSTIA - G16');
colormap(hsub, parula);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);


hsub = subplot(3,3,6);
pcolor(LON, LAT, var_diff.OSTIAfromGPB); shading flat;
shading flat; 
%axis([-59, -48, 8 18]);
caxis([0 0.055]);
title('variance: GPBLEND - OSTIA');
colormap(hsub, parula);
colorbar(hsub);
axis('equal');
axis([-59, -48, 8 18]);


% plot time series:
subplot(3,3,7);
yyaxis left
plot(mean_diff_ts.GPBfromG16,'-','marker','*');
ylabel('spatial diff.');
ylim([-0.2 , 0.2]);
yyaxis right
plot(var_diff_ts.GPBfromG16,'-','marker','*');
ylabel('spatial variance');
ylim([0.005 0.055]);

subplot(3,3,8);
yyaxis left
plot(mean_diff_ts.OSTIAfromG16,'-','marker','*');
ylabel('spatial diff.');
ylim([-0.2 , 0.2]);

yyaxis right
plot(var_diff_ts.OSTIAfromG16,'-','marker','*');
ylabel('spatial variance');
ylim([0.005 0.055]);

subplot(3,3,9);
yyaxis left
plot(mean_diff_ts.OSTIAfromGPB,'-','marker','*');
ylabel('spatial diff.');
ylim([-0.05 0.05]);
yyaxis right
plot(var_diff_ts.OSTIAfromGPB,'-','marker','*');
ylabel('spatial variance');
ylim([0 0.015]);

xc_savefig(gcf,'Figs/compare_L4products', ['movmean_' num2str(wz) 'd_comparison.jpg'],[0 0 12 6]);


%% how variable is the hourly SST from the daily mean GPBlend-L4 sst?

% decide to look at the gpblend L4 product from the brief analysis above:
matFN3 = '_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N';
gpbdata3 = [dataroot filesep 'gpblend' filesep 'gpblend' matFN3];
load(gpbdata3);
gpb = data;
clear data

% apply the blob search function to the L4 data; ()
% what is the large scale mean: (I do need a definition of that).

% try daily mean, large scale>600km:
SSTgap_thres = 0.5;         % when cloudy pixel is less than 50% of the scene, fill the SST map;
dT_thres=0.1;               % temperature anomaly threshould.
cldcv_thres = 0.1;          % threshold to constrain whether or not a warm blob identified is under strong influence of cloud overhead (a.k.a., spatial interpolation.).
minArea=25;                 % 5x5 pixels: 10km (1pixel is 2km across)
dT_thres_arry=[0:0.1:0.8];        % use different threshold to get all the blobs out.

cutoff_scale_km = 1000;

for it = 1:length(gpb.time)
    SST0 = gpb.analysed_sst(:,:,it)';
    lon = double(gpb.lon(:,it));
    lat = double(gpb.lat(:,it));
    
    if mod(it, 10)==0
        checkflag = true;
    else
        checkflag = false;
    end
    SST_LPF(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);    % the size is different here.
    
    SST_anom_spatial = SST0 - SST_LPF(:,:,it);
    fake_cloud_flag = isnan(SST_anom_spatial);
    cloud_flag = fake_cloud_flag;
    
    [LON, LAT] = meshgrid(lon, lat);
    [blob_info, exitflag] = find_SST_blobs(SST_anom_spatial, dT_thres, minArea, cloud_flag, 'cloudcoverage_thres',cldcv_thres, ...
        'LON_grid', LON, 'LAT_grid', LAT,'checkflag',checkflag);
    caxis([-0.8, 0.8]);
    colormap(redblue);
    
    if checkflag
        pause
    end
    
end

% look at the results.
% looked promising! (I can use that to find the cloud mask for now.) % need
% surface wind speed as well from ERA dataset, to-be extracted from PSL.);
