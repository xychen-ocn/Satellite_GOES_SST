% investigating the probability of cloudiness during the ATOMIC period.
% (Can potentially extend to the entire winter seasons (DJF), and for
% multiple years.
% - check Jan and Feb first, use ERA5 to get the mean surface wind speed
% for these period. group days by surface wind.
% check to see if cloudiness freq. is a function of surface wind.
% use the surface wind bins from Bony et al. 2020 paper.

clear all; clc; 
addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');


%% redownload data for the entire Jan and Feb;
% datasvdir = pwd;
% days_requested = 1:59;  % first 2 months.
% download_GOES_SSTdata(days_requested, datasvdir);   % only download data that doesn't exist in the datasvdir.

%% load in all data:
% code design: 

% load in all images:
% datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/g16';
% svdataFN = 'ATOMIC_DFJ_winterseasons2019_2020_SSTmaps.mat';

datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
svdataFN = 'ATOMIC_JanFeb_SSTmaps.mat';


blob_PDF_file = 'probability_density_function_for_warmblob_detection_ATOMIC.mat';
pdf_data = load([datadir filesep blob_PDF_file]);

% load the warm anomaly sampled by RHB and to plot it over the gravel
% day-time cloudiness map.
load([datadir filesep 'blob_data/g16_5d-movingMean_LPF600km/Jan09_RHB_sampled_blobs.mat']);

if ~exist([datadir filesep svdataFN],'file')
    
    files = dir([datadir filesep 'g16/GOES_SST_0*.mat']);
    filenames = {files.name};
    
     % stack up the SST data:
    SST_all = [];
    time_all = [];
    for i = 1:length(filenames)
        matFN = filenames{i};
        tmp = strsplit(matFN,'_');
        dateID = tmp{4}(1:5);
        disp(['getting ' dateID])
        load([datadir filesep 'g16/' matFN]);
        
       
        if i==1
            it0 = 1;
        end
        nt = size(GOES_ATOMIC.sea_surface_temperature,3);
        itN = it0+nt-1;
        SST_all(:,:,it0:itN) = GOES_ATOMIC.sea_surface_temperature;
        time_all(it0:itN) = GOES_ATOMIC.time_num;
        it0 = itN + 1;
    end
    
    lon=GOES_ATOMIC.lon(:,1); lat = GOES_ATOMIC.lat(:,1);
    [LON ,LAT]= meshgrid(lon, lat);
    SST_all = permute(SST_all, [2, 1,3]);
    
    %% QC data:
    % through away maps that is 100% garbage;
    idx_keep = [];
    for i= 1:size(SST_all,3)
        SST_tmp = SST_all(:,:,i);
        idxs = find(isnan(SST_tmp));
        if length(idxs) == numel(SST_tmp)
%             idx_keep = [idx_keep, i];
%         else
            SST_all(:,:,i) = [];
            time_all(:,:,i) = [];
        end
    end
    
%     SST_all_new = SST_all(:,:, idx_keep);
%     SST_all = SST_all_new;
%     
%     time_all = time_all(idx_keep);
    
    nSSTmaps = size(SST_all,3);

    save([datadir filesep svdataFN],'SST_all','time_all','lon','lat','-v7.3');

    
else
    load([datadir filesep svdataFN]);
end

 [SST_LS_trend_spec,~] = estimate_ATOMIC_period_mean_SST_LSG('method','spectrum', 'CutoffScale', 1000);
 %[SST_LS_trend_fit,~] = estimate_ATOMIC_period_mean_SST_LSG;

 % compute gradient:
%  [SST_LSGX_spec, SST_LSGY_spec] = gradient(SST_LS_trend_spec);
%  [SST_LSGX_fit, SST_LSGY_fit] = gradient(SST_LS_trend_fit);
%  LSG_mag_spec = sqrt(SST_LSGX_spec.^2 + SST_LSGY_spec.^2);
%  LSG_mag_fit = sqrt(SST_LSGX_fit.^2 + SST_LSGY_fit.^2);

lon=GOES_ATOMIC.lon(:,1); lat = GOES_ATOMIC.lat(:,1);
[LON ,LAT]= meshgrid(lon, lat);

atomic_region = [-59, 46; 8, 20];
lonmask = lon>=atomic_region(1,1) & lon<=atomic_region(1,2);
latmask = lat>=atomic_region(2,1) & lat<=atomic_region(2,2);




% sort by time:
[time_all, sid] = sort(time_all,'ascend');
SST_all = SST_all(:,:,sid);
       
JanMask = time_all<datenum(2020,2,1) & time_all>=datenum(2020,1,1);
FebMask = time_all>=datenum(2020,2,1);
DecMask = time_all<datenum(2020,1,1);
   
% compute cloud frequency for Jan and Feb, and find the average:
cloudy_flag_all = isnan(SST_all);
JanCount = length(find(JanMask==1));
FebCount = length(find(FebMask==1));
DecCount = nSSTmaps - JanCount-FebCount;

JanCF = sum(double(cloudy_flag_all(:,:,JanMask)),3)./JanCount;
FebCF = sum(double(cloudy_flag_all(:,:,FebMask)),3)./FebCount;
DecCF = sum(double(cloudy_flag_all(:,:,DecMask)),3)./DecCount;

%% compute CF in a day and do 3-day averaging:
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
CF_3dmovmean_2monthsAve = mean(CF_movmean,3,'omitnan');

figure(1); clf;
parula_sub = parula(14);
pcolor(LON, LAT, CF_3dmovmean_2monthsAve); shading flat; colorbar;
hold on
colormap(parula_sub);
[C,F]=contour(LON, LAT, CF_3dmovmean_2monthsAve,[0.3:0.1:0.5],'k');
clabel(C,F,[0.3:0.1:0.5],'color','k');
caxis([0.1 0.8]);
xlim([-59 -48]);
ylim([8 18]);
xlabel('longitude');
ylabel('latitude');
set(gca,'fontsize',14);
title('ave of the 3-day running mean of daily CF');


%%
% load in another data file to condition the map by surface wind speed and
% LTS:
eraFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
era.time = ncread(eraFN,'time');
era.timenum = datenum('1800-01-01') + era.time/24;
era.sugar_tag = ncread(eraFN,'sugar_tag');
era.fish_tag = ncread(eraFN,'fish_tag');
era.gravel_tag = ncread(eraFN,'gravel_tag');
era.flowers_tag = ncread(eraFN,'flowers_tag');

%% compute cloudiness freq from all the 2037 DJF SST maps:
cloudy_flag_all = isnan(SST_all);
% now, for each grid point, check how many times the pixel is cloudy;
cloudy_counts = sum(double(cloudy_flag_all), 3);
%cloudy_freq =cloudy_counts./nt_total;
cloudy_freq_all = cloudy_counts./nSSTmaps;        % 

subplot(1,2,1);
pcolor(LON,LAT, cloudy_freq_all');shading flat; colorbar; title('DJF')
subplot(1,2,2);
pcolor(LON, LAT, ( JanCF+FebCF+DecCF)'./3); shading flat; colorbar; title('mean of DJF');    %this is pretty close to the CF computed from the entire period.


figure
pcolor(LON,LAT, cloudy_freq_all - 0.5*( JanCF+FebCF)); shading flat; colorbar; title('all - Jan-Feb mean');
caxis([-0.1, 0.1]);


% select a type of cloud condition for averaging map:
CF_ref = CF_3dmovmean_2monthsAve;
cloud_types = {'gravel','flowers','sugar','fish'};
for tt = 1:4
    CloudType = cloud_types{tt};
    CN = [CloudType, '_tag'];
    sel_time = era.timenum(era.(CN)==1);
    
    idx_sel=[];
    for it = 1:length(sel_time)
        disp(['era at ', datestr(sel_time(it))]);
        % if the satellite time (hrly) is within the 3hrly window of the
        % era selected time:
        %idx0 = find(time_stack==sel_time(it));
        idxs = find((time_all-sel_time(it)>=0) & (time_all-sel_time(it)<3/24));
        if length(idxs)~=0
            disp(['selected satellite at ', datestr(time_all(idxs(1)))])
            disp(['to satellite at ', datestr(time_all(idxs(end)))])
            
            idx_sel = [idx_sel, idxs];
        end
    end
    idx_sel = unique(idx_sel);  % just in case repeated indices.
    
    grouped_data.SST.(CloudType) = SST_all(:,:,idx_sel);
    grouped_data.time.(CloudType) = time_all(idx_sel);
    
    
    cloudy_flag = isnan(SST_all(:,:,idx_sel));
    nt = length(idx_sel);
    
    % get a sense of the percentage of cloudy pixels for each condition; (which could tell us that perhaps certain condition over predict?
    cloud_spatial_coverage(tt) = mean(sum(double(cloudy_flag), [1,2])./numel(LON));
    
    % now, for each grid point, check how many times the pixel is cloudy;
    cloudy_counts = sum(double(cloudy_flag), 3);
    %cloudy_freq =cloudy_counts./nt_total;
    cloudy_freq.(CloudType) = cloudy_counts./nt;
    cloudy_freq_anom.(CloudType) = (cloudy_counts./nt - CF_ref) ./CF_ref;
    sample_info.sampleSize(tt) = nt;
    sample_info.samplePrc(tt) = nt/nSSTmaps;
    
    
    
    %% further group data into daytime and nighttime samples: (can be done for each type)
    % input will be a collection of SST maps:
    [night.(CloudType), day.(CloudType)]=group_data_into_daytime_and_nighttime_samples(SST_all(:,:,idx_sel), time_all(idx_sel));
    night.(CloudType) = compute_cloudfreq_and_anom(night.(CloudType).SSTmaps, CF_ref);
    day.(CloudType) = compute_cloudfreq_and_anom(day.(CloudType).SSTmaps, CF_ref);
    
    % compare the daytime, nighttime CF anomalies.
end
save([datadir filesep 'cloudiness_freq_ATOMIC.mat'],'cloudy_freq','cloudy_freq_anom','sample_info','cloudy_freq_all','LON','LAT','cloud_types','day','night');


%     freq_lev = [0:0.1:1.1];
%     parula_mask = parula(length(freq_lev)-1);
%     parula_mask(end,:) = [0.2,0.2, 0.2];
% use a warm-blue map:
turbo_all = turbo(64);
parula_all = parula(64);
parula_mid = parula_all(16-8:48+8,:);
turbo_comb = vertcat(turbo_all(24:-1:1,:), parula_all(1:1:48+16,:), turbo_all(64-26:end,:));
% temporary exploration:
for tt = 1:4    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    prc = sample_info.samplePrc(tt)*100;
    
    figure(tt);
   % pcolor(LON, LAT,cloudy_freq.(CloudType)); shading flat;
    pcolor(LON, LAT,cloudy_freq_anom.(CloudType)); shading flat;
    %colormap((parula_mask));
    if tt==3
        colormap(turbo_comb);
    else
    colormap(parula_mid);
    end
    %colormap(redblue(20));
    hb=colorbar;
    %hb.Ticks=[0.1:0.1:0.8];
    if tt~=3
        caxis([-0.5, 0.5]);
    else
        caxis([-1.5, 1.5]);
    end
    set(get(hb,'xlabel'),'string','cloudiness frequency','fontsize',12);
    %set(hb,'Ticks',[0:0.1:1]);
    hold on;
    [C,h]=contour(LON,LAT,cloudy_freq_anom.(CloudType), [-2:0.2:2],'k');
     clabel(C,h, [-2:0.2:2],'labelspacing',800,'color','k','fontsize',12); %,'fontweight','bold');
%     if tt ==4
%         [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
%         clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
%     end
%     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
%     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
    xlabel('longitude');
    ylabel('latitude');
    hold off;
    set(gca,'fontsize',14);
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    title({'Cloudiness Frequency'; [CloudType ' (' ...
        num2str(nt) ';', num2str(prc,'%3.0f')  '%)']},'fontsize',16);
    axis('equal');
   % xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_' CloudType '_v1.jpg'],[0 0 10 8]);
    %xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_'
    %CloudType '_v2.jpg'],[0 0 10 8]);  with SST:

end

% plot the cloud frequency anomaly from the Jan-Feb mean CF in one map:
pos = customize_subplot_size(2,2,0.15, 0.1);
figure(1); clf;
for tt = 1 :4
    % plot anomaly field:
    hsub(tt) = subplot(2,2,tt);
    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    prc = sample_info.samplePrc(tt)*100;
    
    hm = pcolor(LON, LAT,cloudy_freq_anom.(CloudType)); shading flat;
    %colormap((parula_mask));
    if strcmp(CloudType, 'flowers')
        colormap(hsub(tt), turbo_comb);
        caxis([-1.5, 1.5]);
    else
        colormap(hsub(tt), parula_mid);
        caxis([-0.5, 0.5]);
    end
    %colormap(redblue(20));
    hb=colorbar;
    set(get(hb,'xlabel'),'string','anomaly fraction','fontsize',14);
    if ~strcmp(CloudType, 'flowers')
        set(hb,'ticks',[-0.5:0.1:0.5]);
    end
    hold on;
    [C,h]=contour(LON,LAT,cloudy_freq_anom.(CloudType), [-2:0.2:2],'k');
     clabel(C,h, [-2:0.2:2],'labelspacing',800,'color','k','fontsize',10); %,'fontweight','bold');
%     if tt ==4
%          [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
%          clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
%     end
%     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
%     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
    xlabel('Longitude (^{\circ}E)');
    ylabel('Latitude (^{\circ}N)');
    hold off;
    set(gca,'fontsize',14, 'ytick',[8:2:18], 'TickDir','both');
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    title({[CloudType ' (' ...
        num2str(nt) ';', num2str(prc,'%3.0f')  '%)']},'fontsize',16);
    axis('equal');
    set(gca,'pos',pos{tt});
end
xc_savefig(gcf,'Figs',['Anomaly_cloudiness_frequency_in_4Conditions_relative_to_JanFeb_meanstate.jpg'],[0 0 14 10]);


% check differences between night time and day time CF anomalies:
RHB_blobs_time = [RHB_sampled_blobs.time];
ids_night = find(RHB_blobs_time >=datenum('2020-01-08 22:00','yyyy-mm-dd HH:MM'));
ids_day = find(RHB_blobs_time <=datenum('2020-01-09
[~,id] = min(abs(RHB_blobs_time - datenum('2020-01-09 15:00','yyyy-mm-dd HH:MM')));

for tt = 1 %:4
    figure(tt); clf;
    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    %prc = sample_info.samplePrc(tt)*100;
    
    
    for k = 1:2
        if k==1
            val = day.(CloudType).CF_anom;
            titlestr = [CloudType ' condition: daytime'];
                    nn = size(day.(CloudType).SSTmaps,3);

            
        else
            val = night.(CloudType).CF_anom;
            titlestr = [CloudType ' condition: nighttime'];
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
            colormap( parula_mid);
            caxis([-0.5, 0.5]);
        end
        %colormap(redblue(20));
        hb=colorbar;
        set(get(hb,'xlabel'),'string',{'anomaly fraction'},'fontsize',14);
        if ~strcmp(CloudType, 'flowers')
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
        title({titlestr; ['(' num2str(nn) ';', num2str(prc,'%3.0f')  '%)']},'fontsize',16);
        axis('equal')
        ylim([8 18])
        xlim([-59 -48])
        %set(gca,'pos',pos{tt});
    end
    
    %xc_savefig(gcf,'Figs',['Anomaly_cloudiness_frequency_in_day_night_comparison_' CloudType '_addedRHB_blobs.jpg'],[0 0 14 7]);

end


% These four figures are a bit interesting.  (write a bit about this
% tomorrow to describe each figure.
% need to see whether or not any of them are aligned with some eddy
% location.

%% 2nd version: Fish add SST contours on it. (done)



%% 3nd version: Add probaility map during the entire period.
% pdf_levs{1} = [0.005:0.0025:0.01];
% pdf_levs{2} = [0.01:0.0025:0.02];
% pdf_levs{3} = [0.02:0.0025:0.03];
% pdf_color0 = jet(64);

pdf_levs{1} = [0.005:0.0025:0.0075];
pdf_levs{2} = [0.0075:0.0025:0.01];
pdf_levs{3} = [0.01:0.0025:0.0125];
pdf_levs{4} = [0.0125:0.0025:0.015];
pdf_levs{5} = [0.015:0.0025:0.0175];
pdf_levs{6} = [0.0175:0.0025:0.02];
pdf_levs{7} = [0.02:0.0025:0.0225];
pdf_levs{8} = [0.0225:0.0025:0.025];

graysub = gray(64);
pdf_color0 = jet(64);
cindx = floor(linspace(5,50, length(pdf_levs)));
pdf_color = pdf_color0(cindx,:);


%figure(10);
for tt = 1:4
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    prc = sample_info.samplePrc(tt)*100;
    
%     if tt~=3
%         cindx = [1, 28, 40];
%         pdf_color = pdf_color0(cindx,:);
%         
%     else
%         cindx = [1, 40, 62];
%         pdf_color = pdf_color0(cindx,:);
%         
%     end
%     
    figure(tt);clf;
    % subplot(2,2,tt);
    pcolor(LON, LAT,cloudy_freq.(CloudType)); shading flat;
    %colormap((parula_mask));
    colormap(graysub(floor(linspace(20,64,14)),:));
    hb=colorbar;
    %hb.Ticks=[0.1:0.1:0.8];
    caxis([0.1, 0.8]);
    set(get(hb,'xlabel'),'string','cloudiness frequency','fontsize',12);
    %set(hb,'Ticks',[0:0.1:1]);
    hold on;
    [C,h]=contour(LON,LAT,cloudy_freq.(CloudType), [0.1:0.1:0.8],'color',[0.25 0.25 0.25]);
    clabel(C,h, [0.1:0.2:0.8],'labelspacing',800,'color','k','fontsize',12);
    
    for i = 1:length(pdf_levs)
        c = pdf_color(i,:);
        % plot with the conditioned pdf:
        [Cpdf, hpdf] = contour(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_conditioned.(CloudType),pdf_levs{i},'color',c,'linewidth',1.1, 'linestyle','--');
        % now, plot it with the total pdf;
        %[Cpdf, hpdf] = contour(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_all',pdf_levs{i},'color',c,'linewidth',1.1, 'linestyle','--');

        clabel(Cpdf,hpdf,'color',c,'labelSpacing',400);
        hold on;
    end
    %     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
    %     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
    xlabel('longitude');
    ylabel('latitude');
    hold off;
    set(gca,'fontsize',14);
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    title({'Cloudiness Frequency'; [CloudType ' (' ...
        num2str(nt) '; ', num2str(prc,'%3.0f')  '%)']},'fontsize',16);
    axis('equal');
    %xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_'
    %CloudType '_v2.jpg'],[0 0 10 8]);  with SST:
    xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_' CloudType '_with_warmblob_pdf_v2_graymap.jpg'],[0 0 10 8]);
    %xc_savefig(gcf,'Figs',[CloudType '_cloudiness_frequency_GOES_SSTretrival_JanFeb2020_with_warmblob_pdf_alltime.jpg'],[0 0 10 8]);
    
end

pos = customize_subplot_size(2,3,0.15,0.08);

% plot the pdf as shaing as well for reference.
figure(10);clf;
pids = {1, 2, 4, 5, 6};
for tt = 1:5
    
    hsub =subplot(2,3,pids{tt});
    if tt <5
        CloudType = cloud_types{tt};
        
        pcolor(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_conditioned.(CloudType));
    else
        pcolor(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_all');
    end
    shading flat;
    colormap(gray(12));
    hb = colorbar;
    set(get(hb,'xlabel'),'String','Probability Density Estimate');
    caxis([0 0.03]);
    xlabel('longitude');
    ylabel('latitude');
    hold off;
    set(gca,'fontsize',14);
    set(gca,'TickDir','both');
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    if tt <5
        title(CloudType,'fontsize',14);
    else
        title('All', 'fontsize',14);
    end
    xlim([-59, -48]);
    ylim([8, 18]);
    axis('square');
   % if tt <5
        set(hsub, 'pos',pos{pids{tt}});
    %end
    
end
xc_savefig(gcf,'Figs',['Warmblob_pdf_conditioned_by_cloudType_allInOne.jpg'],[0 0 12 8]);
    
    
%% plot the total cloudiness frequency without distinguishing types: add pdf contours:
cloudy_flag_all = isnan(SST_all);

% now, for each grid point, check how many times the pixel is cloudy;
cloudy_counts = sum(double(cloudy_flag_all), 3);
%cloudy_freq =cloudy_counts./nt_total;
cloudy_freq_all = cloudy_counts./nSSTmaps;

%pdf_levs = [0.005:0.0025:0.03];
pdf_levs{1} = [0.005:0.0025:0.0075];
pdf_levs{2} = [0.0075:0.0025:0.01];
pdf_levs{3} = [0.01:0.0025:0.0125];
pdf_levs{4} = [0.0125:0.0025:0.015];
pdf_levs{5} = [0.015:0.0025:0.0175];
pdf_levs{6} = [0.0175:0.0025:0.02];
pdf_levs{7} = [0.02:0.0025:0.0225];
pdf_levs{8} = [0.0225:0.0025:0.025];

graysub = gray(64);
pdf_color0 = jet(64);
cindx = floor(linspace(5,50, length(pdf_levs)));
pdf_color = pdf_color0(cindx,:);

figure;
pcolor(LON, LAT,cloudy_freq_all); shading flat;
%colormap((parula_mask));
colormap(graysub(floor(linspace(20,64,14)),:));
%colormap(gray(14));
hb=colorbar;
%hb.Ticks=[0.1:0.1:0.8];
caxis([0.1, 0.8]);
set(get(hb,'xlabel'),'string','cloudiness frequency','fontsize',12);
%set(hb,'Ticks',[0:0.1:1]);
hold on;
[C,h]=contour(LON,LAT,cloudy_freq_all, [0.1:0.1:0.8],'color',[0.25,0.25,0.25]);
clabel(C,h, [0.1:0.2:0.8],'labelspacing',800,'color','k','fontsize',12);
for i = 1:length(pdf_levs)
   % cindx = [1, 28, 40];
    
    c = pdf_color(i,:);
    [Cpdf, hpdf] = contour(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_all',pdf_levs{i},'color',c,'linewidth',1.1, 'linestyle','--');    
    clabel(Cpdf,hpdf,'color',c,'labelSpacing',400);
end
hold on;
xlabel('longitude');
ylabel('latitude');
hold off;
set(gca,'fontsize',14);
title({'Cloudiness Frequency'; ['All (', num2str(nSSTmaps)  ')']},'fontsize',16);
axis('equal');
xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_ALL_with_warmblob_probability_v3.jpg'],[0 0 10 8]);




%% I am a bit curious if I randomly select the same amount of data for each group, will I get localized increase in CF? 
%% if so? then maybe what I have been showing is nonsense?

% randomly select samples and compare to the meanstate.
%ridx = [];

for it = 1:10
    SST_all_operate = SST_all;
    time_all_operate = time_all;
    N0 = 1375;
    for i = 1:4
        if i==1
            N = N0;
        end
        ridx = randperm(N, sample_info.sampleSize(i));
        CF_type(i,it).SST_maps=SST_all_operate(:,:,ridx);
        CF_type(i,it).time = time_all_operate(ridx);
        CF_type(i,it).idx  = ridx;
        SST_all_operate(:,:,ridx)=[];       % remove the sample already drew.
        time_all_operate(ridx)=[];
        N = size(SST_all_operate,3);
        
        % compute cloudiness frequency:
        cloud_flag = isnan(CF_type(i,it).SST_maps);
        CF_type(i,it).CF = sum(double(cloud_flag),3)./size(cloud_flag,3);
        CF_type(i,it).CFa = (CF_type(i,it).CF  - cloudy_freq_all)./cloudy_freq_all;
        
    end
end


figure(11);
for i = 1:4
    subplot(2,2,i);
    for it =1 :10
        CF_type_all(:,:,it) = CF_type(i,it).CFa;
    end
    ave_CFa = mean(CF_type_all,3);
    pcolor(LON, LAT, ave_CFa);
    shading flat;
    hold on;
    [C,h]=contour(LON,LAT,ave_CFa, [-2:0.2:2],'k');
    clabel(C,h,[-0.8:0.2:0.8],'color','k');
    
    colorbar;
    caxis([-0.5  0.5]);
    title({['random group ' num2str(i)]; [num2str(sample_info.sampleSize(i))]});

end
xc_savefig(gcf,'Figs','check_cloudiness_freq_changes_in_randomly_separated_groups_ave_of_10assignments.jpg',[0 0 12, 10]);


for it = 1:10
    figure(12);clf
    
    for i = 1:4
        subplot(2,2,i);
        
        pcolor(LON, LAT, CF_type(i,it).CFa);
        shading flat;
        hold on;
        [C,h]=contour(LON,LAT,CF_type(i,it).CFa, [-2:0.2:2],'k');
        clabel(C,h,[-0.8:0.2:0.8],'color','k');
        
        colorbar;
        caxis([-0.5  0.5]);
        title({['random group ' num2str(i)]; [num2str(sample_info.sampleSize(i))]});
        
        
    end
    xc_savefig(gcf,'Figs',['check_cloudiness_freq_changes_in_randomly_separated_groups_assignmentNum' num2str(it,'%2.2i') '.jpg'],[0 0 12, 10]);
    
    
end