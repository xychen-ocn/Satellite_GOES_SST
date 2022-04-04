% This script is used to compute cloudiness frequency differently;
%
clear all; clc; close all;

datadir = pwd
load('ATOMIC_JanFeb_SSTmaps.mat');

% load ncfile;
eraFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
era.time = ncread(eraFN,'time');
era.timenum = datenum('1800-01-01') + era.time/24;
era.U10_domainMean = ncread(eraFN,'U10_domainAve');
era.LTS_domainMean = ncread(eraFN, 'LTS_domainAve');

% do daily average on two quantity
eradays = unique(floor(era.timenum));
for i = 1:length(eradays)
    tmask = era.timenum>=eradays(i) & era.timenum<eradays(i)+1;
    era.U10_dailyDomainAve(i) = mean(era.U10_domainMean(tmask));
    era.LTS_dailyDomainAve(i) = mean(era.LTS_domainMean(tmask));
end
era.day_timenum = eradays;

figure(1); clf;
%  3hrly data:
h(2)= plot(era.LTS_domainMean, era.U10_domainMean, '.','linestyle','none','color',[0.5 0.5 0.5]);
hold on;
h(1)= plot(era.LTS_dailyDomainAve, era.U10_dailyDomainAve,'.b','markersize',18);

plot([15, 15], [4, 12],'--r','linewidth',1.2);
plot([11, 19], [8, 8],'--r','linewidth',1.2);
%title('Large Scale Atmospheric conditions in Jan, Feb 2020')
set(gca,'fontsize',16);
xlabel('LTS (K)');
ylabel('U_{10} (m/s)');
grid on
xlim([11, 19]); ylim([4 12]);
legend(h, {'daily averaged','3-hourly'});
xc_savefig(gcf,'./','ERA5_U10_LTS_dailyave_JanFeb2020.jpg',[0 0 8 6]);

% use the daily average quantity to provide a condition tag;
cond_U10 = era.U10_dailyDomainAve>=8;
cond_LTS = era.LTS_dailyDomainAve>=15;

cond.flowers = cond_U10& cond_LTS;
cond.gravel = cond_U10 & ~cond_LTS;
cond.sugar = ~cond_U10 & ~cond_LTS;
cond.fish = ~cond_U10 & cond_LTS;

[LON, LAT]=meshgrid(lon,lat);

cloud_types = {'gravel','flowers','sugar','fish'};

turbo_all = turbo(64);
parula_all = parula(64);
parula_mid = parula_all(16-8:48+8,:);
turbo_comb = vertcat(turbo_all(24:-1:1,:), parula_all(1:1:48+16,:), turbo_all(64-26:end,:));


SSTmapTypetag = zeros(size(time_all));
for m= 1:4
    CT = cloud_types{m};
    dayType.(CT) = era.day_timenum(cond.(CT));
    
    % now: find the 3d and 5d running mean for each day type of large-scale
% condition.
    % find SST maps within the selected days
    idx_sel = [];
    for id = 1:length(dayType.(CT))
        dhere = dayType.(CT)(id);
        tmask_1d = time_all>=dhere & time_all<dhere+1;
        idx_sel = [idx_sel, find(tmask_1d==1)];

    end
    
    SSTmapTypetag(idx_sel)=m;
end

for m= 1:4
    CT = cloud_types{m};
    dayType.(CT) = era.day_timenum(cond.(CT));
    
    % now: find the 3d and 5d running mean for each day type of large-scale
% condition.
    % find SST maps within the selected days
    idx_sel = [];
    for id = 1:length(dayType.(CT))
        dhere = dayType.(CT)(id);
        tmask_1d = time_all>=dhere & time_all<dhere+1;
        tmask_3d = time_all>=dhere-1 & time_all<dhere+2;
        tmask_5d = time_all>=dhere-2 & time_all<dhere+3;
        
        idx_sel = [idx_sel, find(tmask_1d==1)];
        
        SSTmaps_1d = SST_all(:,:,tmask_1d);
        SSTmaps_3d = SST_all(:,:,tmask_3d);
        SSTmaps_5d = SST_all(:,:,tmask_5d);
        
        nt_1d(id) = size(SSTmaps_1d,3);
        nt_3d(id) = size(SSTmaps_1d,3);
        nt_5d(id) = size(SSTmaps_1d,3);

        if nt_1d(id)>1
            daily_CF(:,:,id) = compute_cloudfreq(SSTmaps_1d);
        else
            daily_CF(:,:,id)=NaN;
        end
        CF_3day(:,:,id) = compute_cloudfreq(SSTmaps_3d);
        CF_5day(:,:,id) = compute_cloudfreq(SSTmaps_5d);

    end
    
    SSTmapTypetag(idx_sel)=m;
    
    
    CF_3danom = (daily_CF - CF_3day)./CF_3day;
    CF_5danom = (daily_CF - CF_5day)./CF_5day;

    CF.(CT).daily_CF_ensm = mean(daily_CF,3,'omitnan');
    CF.(CT).CF_3danom_ensm = mean(CF_3danom,3,'omitnan');
    CF.(CT).CF_5danom_ensm = mean(CF_5danom,3,'omitnan');
    
    mean_CF3d = mean(CF_3day,3,'omitnan');
    mean_CF5d = mean(CF_5day,3,'omitnan');
    
    parula_all = parula(64);
    parula_mid = parula_all(16-10:48+10,:);
    
    figure(2+k)
    subplot(2,1,1)
    pcolor(lon, lat, CF.(CT).CF_3danom_ensm);
    shading flat; colorbar;
    colormap(parula_mid); hold on;
    [C,h]=contour(lon, lat, CF.(CT).CF_3danom_ensm,[-0.5:0.1:0.5],'k');
    caxis([-0.5 0.5]);
    xlim([-59 -48]);
    ylim([8 18]);
    title({['averaged CF anomalies in ' CT '-favored condition']; 'relative to  3-day running mean'})
    
    subplot(2,1,2)
    pcolor(lon, lat, CF.(CT).CF_5danom_ensm);
    shading flat; colorbar; hold on;
    colormap(parula_mid);
    [C,h]=contour(lon, lat, CF.(CT).CF_3danom_ensm,[-0.5:0.1:0.5],'k');
    caxis([-0.5 0.5]);
    xlim([-59 -48]);
    ylim([8 18]);
    title({['averaged CF anomalies in ' CT '-favored condition']; 'relative to  5-day running mean'})

    pause;
    
    
    % further separate into daytime and nighttime anomaly:
    %% further group data into daytime and nighttime samples: (can be done for each type)
    % input will be a collection of SST maps:
    [night.(CT), day.(CT)]=group_data_into_daytime_and_nighttime_samples(SST_all(:,:,idx_sel), time_all(idx_sel));
    night_ref3d.(CT) = compute_cloudfreq_and_anom(night.(CT).SSTmaps, mean(CF_3day,3));
    day_ref3d.(CT) = compute_cloudfreq_and_anom(day.(CT).SSTmaps, mean(CF_3day,3));
    
    night_ref5d.(CT) = compute_cloudfreq_and_anom(night.(CT).SSTmaps,  mean(CF_5day,3));
    day_ref5d.(CT) = compute_cloudfreq_and_anom(day.(CT).SSTmaps,  mean(CF_5day,3));
    
    
    % plot day and night:
    for j= 1:2
        nt = length(idx_sel);
        if j==1
            day_data = day_ref3d;
            night_data = night_ref3d;
        else
            day_data = day_ref5d;
            night_data = night_ref5d;
        end
        
        for k = 1:2
            if k==1
                val = day_data.(CT).CF_anom;
                titlestr = [CT ' condition: daytime'];
                nn = size(day.(CT).SSTmaps,3);
                
                
            else
                val = night_data.(CT).CF_anom;
                titlestr = [CT ' condition: nighttime'];
                nn = size(night.(CT).SSTmaps,3);
                
            end
            prc = nn/nt*100;
            
            figure(3)
            subplot(2,2,k+2*(j-1));
            
            hm = pcolor(LON, LAT, val); shading flat;
            %colormap((parula_mask));
            if strcmp(CT, 'flowers')
                colormap( turbo_comb);
                caxis([-1.5, 1.5]);
            else
                colormap( parula_mid);
                caxis([-0.5, 0.5]);
            end
            %colormap(redblue(20));
            hb=colorbar;
            set(get(hb,'xlabel'),'string',{'anomaly fraction'},'fontsize',14);
            if ~strcmp(CT, 'flowers')
                set(hb,'ticks',[-0.5:0.1:0.5]);
            end
            hold on;
            [C,h]=contour(LON,LAT,val, [-2:0.2:2],'k');
            clabel(C,h, [-2:0.2:2],'labelspacing',800,'color','k','fontsize',10); %,'fontweight','bold');
            %     if tt ==4
            %         [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
            %         clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
            hold on;
            % plot RHB sampled blobs near 15UTC
            %         for id = 1:5
            %         contour(RHB_sampled_blobs(id).blobImageCoord.lon, RHB_sampled_blobs(id).blobImageCoord.lat, ...
            %                 RHB_sampled_blobs(id).blobImage{1},'m','linewidth',1.1,'linestyle',':');
            %         end
            %     end
            %     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
            %     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
            xlabel('Longitude (^{\circ}E)');
            ylabel('Latitude (^{\circ}N)');
            hold off;
            set(gca,'fontsize',14, 'ytick',[8:2:18], 'TickDir','both');
            %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
            if j ==1
                title({titlestr; ['(' num2str(nn) ';', num2str(prc,'%3.0f')  '%)'];'reference: 3-d running mean'},'fontsize',16);
            else
                title({titlestr; ['(' num2str(nn) ';', num2str(prc,'%3.0f')  '%)'];'reference: 5-d running mean'},'fontsize',16);
                
            end
            axis('equal')
            ylim([8 18])
            xlim([-59 -48])
            %set(gca,'pos',pos{tt});
        end
        
        
    end
    
    
    
    
end


% check mean-3d and mean_5d and compare to the mean from the entire Jan,
% Feb months;
CF_JanFeb = compute_cloudfreq(SST_all);

figure
subplot(1,2,1)
pcolor(LON, LAT, mean_CF3d); shading flat; 
colorbar;
hold on;
colormap(parula(14));
[C,h]=contour(lon, lat, CF.(CT).CF_3danom_ensm,[-0.5:0.1:0.5],'k');
caxis([0.1 0.8]);
xlim([-59 -48]);
ylim([8 18]);
axis('square');
title('composite of 3-day running mean from gravel-condition')
set(gca,'fontsize',14);

subplot(1,2,2)
pcolor(LON, LAT, CF_JanFeb); shading flat;
colorbar;
hold on;
colormap(parula(14));
[C,h]=contour(lon, lat, CF.(CT).CF_3danom_ensm,[-0.5:0.1:0.5],'k');
caxis([0.1 0.8]);
xlim([-59 -48]);
ylim([8 18]);
axis('square');
title('CF from the entire sample set')
set(gca,'fontsize',14);

% compute reference state by moving average the daily CF using a 5-day
% window, and then composite on all the 60 days in Jan, Feb 2020.
data_days = unique(floor(time_all));
    
for id =1:length(data_days)
    
    tmask = time_all>=data_days(id) & time_all<data_days(id)+1;
    cld_flag = isnan(SST_all(:,:,tmask));
    nt = length(find(tmask==1));
    daily_samplesize(id) = nt;
    if nt>1
    CF_daily(:,:,id) = sum(cld_flag,3)./nt;
    else
        CF_daily(:,:,id) = nan(size(LON));
    end
    
end
CF_movmean = movmean(CF_daily,5,3,'omitnan');
CF_5dmovmean_2monthsAve = mean(CF_movmean,3,'omitnan');

%% plot reference state CF with large scale SST overlay on it:
% load in the large scale mean SST use the utility function I wrote:
product='g16'; wz=5; LPF_scale_km=600; search_region = [-59, -45; 8 22];    % avoid fft complication by land during 
SST_mean = get_large_scale_temporal_mean(product, wz, LPF_scale_km, search_region);
% averaged 5-day running mean SST over the Jan, Feb period.
ave_SST_LS = mean(SST_mean.LPF,3);
[SSTLON, SSTLAT] = meshgrid(SST_mean.lon, SST_mean.lat);


figure
%subplot(1,2,1)
pcolor(LON, LAT, CF_5dmovmean_2monthsAve); shading flat; 
hb=colorbar;
hold on;
colormap(parula(14));
[C,h]=contour(lon, lat, CF_5dmovmean_2monthsAve,[0.1:0.1:1],'k');
clabel(C,h,[0.1:0.1:0.8],'color','k')
caxis([0.1 0.8]);
xlim([-59 -48]);
ylim([8 18]);
%% add the large scale SST:
[C3,h3]= contour(SSTLON, SSTLAT, ave_SST_LS-273.15, [25:0.5:28],'w','linestyle','--','linewidth',1.2);
clabel(C3,h3,[25:0.5:28],'color','w', 'labelspacing',300,'fontsize',12);
%axis('square');
hold off
%title('composite of 5-day running mean CF')
set(gca,'fontsize',16);
xlabel('Longitude (^{\circ}E)');
ylabel('Latitude (^{\circ}N)');
set(get(hb,'xlabel'),'string','cloudiness frequency')
xc_savefig(gcf,'Figs','CF_reference_5dayrunningmean.jpg',[0 0 10 8]);


% subplot(1,2,2)
% pcolor(LON, LAT, CF_JanFeb); shading flat;
% colorbar;
% hold on;
% colormap(parula(14));
% [C,h]=contour(lon, lat, CF.(CT).CF_JanFeb,[-0.5:0.1:0.5],'k');
% caxis([0.1 0.8]);
% xlim([-59 -48]);
% ylim([8 18]);
% axis('square');
% %title('Cloudiness Frequency in Jan, Feb')
% set(gca,'fontsize',14);
% 
% 
% 

%% compute cloudiness freq and its change with different cloud conditions:
CF_ref = CF_5dmovmean_2monthsAve;
cloud_types = {'gravel','flowers','sugar','fish'};
for tt = 1:4
    CloudType = cloud_types{tt};
    CN = [CloudType, '_tag'];
    
    
    idx_sel = SSTmapTypetag==tt;  % just in case repeated indices.
    
    cloudy_freq_cond.(CloudType) = compute_cloudfreq_and_anom(SST_all(:,:,idx_sel), CF_ref);

    sample_info.sampleSize(tt) = length(find(idx_sel==1));
    sample_info.samplePrc(tt) = length(idx_sel)/size(SST_all,3);
    
    
    %% further group data into daytime and nighttime samples: (can be done for each type)
    % input will be a collection of SST maps:
    [night.(CloudType), day.(CloudType)]=group_data_into_daytime_and_nighttime_samples(SST_all(:,:,idx_sel), time_all(idx_sel));
    night.(CloudType) = compute_cloudfreq_and_anom(night.(CloudType).SSTmaps, CF_ref);
    day.(CloudType) = compute_cloudfreq_and_anom(day.(CloudType).SSTmaps, CF_ref);
%     
    % compare the daytime, nighttime CF anomalies.
end
save([datadir filesep 'cloudiness_freq_ATOMIC_JanFeb_updated.mat'],'cloudy_freq_cond','sample_info','CF_3dmovmean_2monthsAve','CF_JanFeb','LON','LAT','cloud_types','day','night','-v7.3');

save([datadir filesep 'cloudiness_freq_ATOMIC_JanFeb_updated.mat'],'sample_info','-append');



% check how what the days are:
figure(4); clf;
colorn={'c','r','k','b'};
for m = 1:4
    sel = SSTmapTypetag==m;
    CT = cloud_types{m};
    plot(time_all(sel), SSTmapTypetag(sel),'.','color', colorn{m});
    %plot(dayType.(CT), m, 'x','color',k
    hold on
end
legend(cloud_types);
ylim([0 5]);
set(gca,'ytick',[1:4], 'yticklabel',cloud_types,'xtick',[time_all(1):10:time_all(1)+60])
grid on
datetick('x','mmmdd','keepticks');
%xlim([time_all(1),time_all(1)+59 ])

set(gca,'fontsize',14);






% search for hrly SST maps within a day to compute the CF;
% compare this to a 3day or 5day running mean centered on the day of
% interest;
load([datadir filesep 'cloudiness_freq_ATOMIC_JanFeb_updated.mat']);
% load significance test:
load('sigTest/significant_threshold_from_individual_sampling_in_4regimes.mat');
load('locationsPDF_of_extreme_effective_SSTgrad_more_time_consistent.mat');
siglev = indvtest_extreme_thres(:,2);   % extreme value of the change:
siglev_round = round(siglev*100)/100;


%% plot changes in CF for 4 different cases.
pos = customize_subplot_size(2,2,0.15, 0.1);
%hot_all = hot(64);
pyRdBu = getPyPlot_cMap('RdBu_r');
redblue_all = pyRdBu ;
redblue_mid = pyRdBu(33:2:96,:);
%rb_comb = vertcat(turbo_all(24:-1:1,:), redblue_mid, flipud(turbo_all(64-24+1:end,:)));

figure(2); clf;
for tt = 1 :4
    % plot anomaly field:
    hsub(tt) = subplot(2,2,tt);
    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    Ntot = sum(sample_info.sampleSize);
    prc = nt /Ntot * 100;
    %prc = sample_info.samplePrc(tt)*100;
    
    hm = pcolor(LON, LAT,cloudy_freq_cond.(CloudType).CF_anom); shading flat;
    %colormap((parula_mask));
    if strcmp(CloudType, 'flowers')
        colormap(hsub(tt), pyRdBu);
        caxis([-1, 1]);
    else
        colormap(hsub(tt), redblue_mid);
        caxis([-0.5, 0.5]);
    end
    %colormap(redblue(20));
    hb=colorbar;
    set(get(hb,'xlabel'),'string','fractional change','fontsize',14);
    if ~strcmp(CloudType, 'flowers')
        set(hb,'ticks',[-0.5:0.1:0.5]);
    else
        set(hb,'ticks',[-1:0.2:1]); %, -0.5, 0.5]));
    end
    hold on;
    [C,h]=contour(LON,LAT,cloudy_freq_cond.(CloudType).CF_anom, [0:0.2:2],'color',[0.45 0.45 0.45]);
     %clabel(C,h, [0.2:0.2:2],'labelspacing',800,'color','k','fontsize',10); %,'fontweight','bold');
    [Cn,hn]=contour(LON,LAT,cloudy_freq_cond.(CloudType).CF_anom, [-2:0.2:-0.2],'k','linestyle',':','color',[0.45 0.45 0.45]);
     %clabel(Cn,hn, [-0.8:0.2:-0.2],'labelspacing',800,'color','k','fontsize',10); %,'fontweight','bold');
     
     % add significant level: 
     if strcmp(CloudType, 'flowers')
         [Csig, hsig] = contour(LON, LAT, cloudy_freq_cond.(CloudType).CF_anom, [siglev_round(tt):0.4:1],'-k','linewidth',1.1);
         
         %clabel(Csig, hsig, [siglev_round(tt):0.4:1],'color','k');
         
     else
         [Csig, hsig] = contour(LON, LAT, cloudy_freq_cond.(CloudType).CF_anom, [siglev_round(tt):0.2:1],'-k','linewidth',1.1);
         %clabel(Csig, hsig, [siglev_round(tt):0.2:1],'color','k');
         
     end
     
     %% add large scale SST (March 30, 2022)
     [CSST,hSST]= contour(LON_sub, LAT_sub, mean(SST_LS.(CloudType),3,'omitnan')-273.15, [25:0.5:28],'--b','linewidth',1.2);
     clabel(CSST,hSST, [25:0.5:28]);

     
     %% add probability density estimates for warm/cold fronts in each type (March 30, 2022)
     contour(lon_bin, lat_bin,canom_pdf{tt}, pdf_levs{tt}, 'c','linewidth',1.1);
     contour(lon_bin, lat_bin,wanom_pdf{tt}, pdf_levs{tt}, 'm','linewidth',1.1);

     
     
%     if tt ==4
%          [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
%          clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
%     end
%     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
%     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
    xlabel('Longitude (^{\circ}E)');
    ylabel('Latitude (^{\circ}N)');
    hold off;
    set(gca,'fontsize',14, 'ytick',[8:2:18], 'TickDir','both','TickLength',[0.015,0.01]);
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    title({[CloudType '-favored regime (' ...
        num2str(nt) ';', num2str(round(prc),'%3.0f')  '%)']},'fontsize',16);
    %axis('equal');
    xlim([-59, -48]); ylim([8 18])
    set(gca,'pos',pos{tt});
end
xc_savefig(gcf,'Figs',['Anomaly_cloudiness_frequency_in_4Conditions_relative_to_5daymovingmean_refstate_March31_added_location_of_extreme_SSTgrad.jpg'],[0 0 12 10]);


%% plot changes in frequency for the gravel case; overlay two regions of interests + the locations of warm anomalies sampled by RHB.
% plot track of RHB instead of the blobs it captured, highlight location of
% the peak of the blobs.

RHB_segdir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/object_oriented_scripts/processed_data';
RHB_segFN = 'RHB_along_wind_segments_ready_for_wavelet_coherence_updated_Ro_24-Jan-2022.mat';
load([RHB_segdir filesep RHB_segFN]);
RHB_Jan09seg.lon = RHB_straight_segs(1).along_wind.lon;
RHB_Jan09seg.lat = RHB_straight_segs(1).along_wind.lat;

blobs_time = [datenum('2020-01-09 06:00','yyyy-mm-dd HH:MM'), datenum('2020-01-09 16:30','yyyy-mm-dd HH:MM')];
RHB_Jan09seg.blon = interp1(RHB_straight_segs(1).along_wind.time, RHB_Jan09seg.lon, blobs_time);
RHB_Jan09seg.blat = interp1(RHB_straight_segs(1).along_wind.time, RHB_Jan09seg.lat, blobs_time);

% RHB_blobdataFN = [datadir filesep 'blob_data/g16_5d-movingMean_LPF600km/Jan09_RHB_sampled_blobs.mat'];
% load(RHB_blobdataFN);

% RHB_blobs_time = [RHB_sampled_blobs.time];
% [~, bid(2)] = min(abs(RHB_blobs_time - datenum('2020-01-09 07:00','yyyy-mm-dd HH:MM')));
% [~,bid(1)] = min(abs(RHB_blobs_time - datenum('2020-01-09 18:00','yyyy-mm-dd HH:MM')));

%turbo_all = turbo(64);

pyRdBu = getPyPlot_cMap('RdBu_r');
redblue_all = pyRdBu ;
redblue_mid = pyRdBu(33:2:64+32,:);
%turbo_comb = vertcat(turbo_all(24:-1:1,:), parula_all(1:1:48+16,:), turbo_all(64-26:end,:));

% ROI(1).coverage = [-53, -51; 13.5, 15.5];
% ROI(1).name = 'RHB_sampled_Hotspots';
% 
% ROI(2).coverage = [-56, -54; 10, 12];
% ROI(2).name = 'non_hotspot_area';

Rvx(:,1) = [-55, -55, -51, -51, -55];
Rvy(:,1) = [13.5, 15.5, 15.5, 13.5, 13.5];

Rvx(:,2) = [-56, -56, -54, -54, -56];
Rvy(:,2) = [10, 12, 12, 10, 10];

for tt = 1 %:4
    figure(tt); clf;
    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    %prc = sample_info.samplePrc(tt)*100;
    
    
    for k = 1:2
        if k==1
            val = day.(CloudType).CF_anom;
            titlestr = [CloudType '-favored regime: daytime'];
                    nn = size(day.(CloudType).SSTmaps,3);

            
        else
            val = night.(CloudType).CF_anom;
            titlestr = [CloudType '-favored regime: nighttime'];
                    nn = size(night.(CloudType).SSTmaps,3);

        end
        prc = nn/nt*100;
        
        subplot(1,2,k);
        
        hm = pcolor(LON, LAT, val); shading flat;
        %colormap((parula_mask));
        if strcmp(CloudType, 'flowers')
            colormap( turbo_comb);
            caxis([-1.5, 1.5]);
        else
            colormap(redblue_mid);
            caxis([-0.5, 0.5]);
        end
        %colormap(redblue(20));
        hb=colorbar;
        set(get(hb,'xlabel'),'string',{'fractional change'},'fontsize',14);
        if ~strcmp(CloudType, 'flowers')
            set(hb,'ticks',[-0.5:0.1:0.5]);
        end
        hold on;
        contour(LON,LAT,val, [0:0.2:2],'-','color',[0.5 0.5 0.5]);
        contour(LON,LAT, val,[-2:0.2:-0.2],'--','color',[0.5 0.5 0.5]);
        % highlight significant levels:
        [Csig, hsig] = contour(LON, LAT, val, [0.2:0.2:1],'-k','linewidth',1.2);

        
       % clabel(C,h, 0.4,'labelspacing',800,'color','k','fontsize',10); %,'fontweight','bold');
        %     if tt ==4
%         [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
%         clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
        hold on;
        %% add two regions of interest:
%         if k==1
%             plot(Rvx(:,1), Rvy(:,1),'--','linewidth',2.5,'color','m');
%             text(Rvx(1,1), Rvy(1,1), {'A'}, 'horizontalAlignment','left','fontsize',20,'fontweight','bold');
%         else
%             plot(Rvx(:,1), Rvy(:,1),'--','linewidth',2.5,'color','m');
%             text(Rvx(1,1), Rvy(1,1),'A','horizontalAlignment','left','fontsize',20,'fontweight','bold');
%         end
        
         %% add large scale SST (March 30, 2022)
     [CSST,hSST]= contour(LON_sub, LAT_sub, mean(SST_LS.(CloudType),3,'omitnan')-273.15, [25:0.5:28],'--b','linewidth',1.2);
     clabel(CSST,hSST, [25:0.5:28]);

     
     %% add probability density estimates for warm/cold fronts in each type (March 30, 2022)
     contour(lon_bin, lat_bin,canom_pdf{tt}, pdf_levs{tt}, 'c','linewidth',1.1);
     contour(lon_bin, lat_bin,wanom_pdf{tt}, pdf_levs{tt}, 'm','linewidth',1.1);

        
        % plot RHB sampled blobs near 15UTC
        %% plot RHB trajectory instead.
        plot(RHB_Jan09seg.lon, RHB_Jan09seg.lat,'-.','color','b','linewidth',2);
        % indicate location of the blob:
        plot(RHB_Jan09seg.blon, RHB_Jan09seg.blat,'x','linestyle','none','color','b','linewidth',2,'markersize',10);
        
%         for id = 1:5
%         contour(RHB_sampled_blobs(bid(k)).blobImageCoord.lon, RHB_sampled_blobs(bid(k)).blobImageCoord.lat, ...
%                 RHB_sampled_blobs(bid(k)).blobImage{1},'c','linewidth',1.,'linestyle','-');
% %         end
        %     end
        %     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
        %     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
        xlabel('Longitude (^{\circ}E)');
        ylabel('Latitude (^{\circ}N)');
        hold off;
       % set(gca,'fontsize',14, 'ytick',[8:2:18], 'TickDir','both');
        set(gca,'fontsize',14, 'ytick',[8:2:18], 'TickDir','both','TickLength',[0.015,0.01]);

        %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
        title({titlestr; ['(' num2str(nn) ';', num2str(prc,'%3.0f')  '%)']},'fontsize',16);
        axis('equal')
        ylim([8 18])
        xlim([-59 -48])
        %set(gca,'pos',pos{tt});
    end
    
    xc_savefig(gcf,'Figs',['Anomaly_cloudiness_frequency_in_day_night_comparison_' CloudType '_addedRHB_blobs_v3_added_SSTgradLocs.jpg'],[0 0 12 7]);

end

