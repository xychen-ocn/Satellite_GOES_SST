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
datasvdir = pwd;
days_requested = 1:59;  % first 2 months.
download_GOES_SSTdata(days_requested, datasvdir);   % only download data that doesn't exist in the datasvdir.

%% load in all data:
% code design: 

% load in all images:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
svdataFN = 'ATOMIC_JanFeb_SSTmaps.mat';
blob_PDF_file = 'probability_density_function_for_warmblob_detection_ATOMIC.mat';
pdf_data = load([datadir filesep blob_PDF_file]);

if ~exist([datadir filesep svdataFN],'file')
    
    files = dir([datadir filesep 'GOES_SST*.mat']);
    filenames = {files.name};
    
    
    for i = 1:length(filenames)
        matFN = filenames{i};
        tmp = strsplit(matFN,'_');
        dateID = tmp{end}(1:5);
        disp(['getting ' dateID])
        load([datadir filesep matFN]);
        
        data_all(i)=GOES_ATOMIC;
    end
    
    % stack up the SST data:
    SST_all = [];
    time_all = [];
    for i= 1:length(filenames)
        if i==1
            it0 = 1;
        end
        nt = size(data_all(i).sea_surface_temperature,3);
        itN = it0+nt-1;
        SST_all(:,:,it0:itN) = data_all(i).sea_surface_temperature;
        time_all(it0:itN) = data_all(i).time_num;
        it0 = itN + 1;
    end
else
    load([datadir filesep svdataFN]);
end

 [SST_LS_trend_spec,~] = estimate_ATOMIC_period_mean_SST_LSG('method','spectrum', 'CutoffScale', 1000);
 [SST_LS_trend_fit,~] = estimate_ATOMIC_period_mean_SST_LSG;

 % compute gradient:
%  [SST_LSGX_spec, SST_LSGY_spec] = gradient(SST_LS_trend_spec);
%  [SST_LSGX_fit, SST_LSGY_fit] = gradient(SST_LS_trend_fit);
%  LSG_mag_spec = sqrt(SST_LSGX_spec.^2 + SST_LSGY_spec.^2);
%  LSG_mag_fit = sqrt(SST_LSGX_fit.^2 + SST_LSGY_fit.^2);

%% QC data:
% through away maps that is 100% garbage;
idx_keep = [];
for i= 1:size(SST_all,3)
    SST_tmp = SST_all(:,:,i);
    idxs = find(isnan(SST_tmp));
    if length(idxs) < numel(SST_tmp)
        idx_keep = [idx_keep, i];
    else
        % drop;
    end
end

SST_all_new = SST_all(:,:, idx_keep);
SST_all = SST_all_new;

time_all = time_all(idx_keep);

nSSTmaps = size(SST_all,3);
        
    

% load in another data file to condition the map by surface wind speed and
% LTS:
eraFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
era.time = ncread(eraFN,'time');
era.timenum = datenum('1800-01-01') + era.time/24;
era.sugar_tag = ncread(eraFN,'sugar_tag');
era.fish_tag = ncread(eraFN,'fish_tag');
era.gravel_tag = ncread(eraFN,'gravel_tag');
era.flowers_tag = ncread(eraFN,'flowers_tag');


cloudy_flag_all = isnan(SST_all);

% now, for each grid point, check how many times the pixel is cloudy;
cloudy_counts = sum(double(cloudy_flag_all), 3);
%cloudy_freq =cloudy_counts./nt_total;
cloudy_freq_all = cloudy_counts./nSSTmaps;


% select a type of cloud condition for averaging map:
cloud_types = {'sugar','gravel','flowers','fish'};
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
    
    
    cloudy_flag = isnan(SST_all(:,:,idx_sel));
    nt = length(idx_sel);
    
    % now, for each grid point, check how many times the pixel is cloudy;
    cloudy_counts = sum(double(cloudy_flag), 3);
    %cloudy_freq =cloudy_counts./nt_total;
    cloudy_freq.(CloudType) = cloudy_counts./nt;
    sample_info.sampleSize(tt) = nt;
    sample_info.samplePrc(tt) = nt/nSSTmaps;
end
save([datadir filesep 'cloudiness_freq_ATOMIC.mat'],'cloudy_freq','sample_info','cloudy_freq_all','LON','LAT');
   
%     freq_lev = [0:0.1:1.1];
%     parula_mask = parula(length(freq_lev)-1);
%     parula_mask(end,:) = [0.2,0.2, 0.2];
for tt = 1:4    
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    prc = sample_info.samplePrc(tt)*100;
    
    figure(tt);
    pcolor(LON, LAT,cloudy_freq.(CloudType)); shading flat;
    %colormap((parula_mask));
    colormap(gray(14));
    hb=colorbar;
    %hb.Ticks=[0.1:0.1:0.8];
    caxis([0.1, 0.8]);
    set(get(hb,'xlabel'),'string','cloudiness frequency','fontsize',12);
    %set(hb,'Ticks',[0:0.1:1]);
    hold on;
    [C,h]=contour(LON,LAT,cloudy_freq.(CloudType), [0.1:0.1:0.8],'k');
    clabel(C,h, [0.1:0.1:0.8],'labelspacing',800,'color','k','fontsize',12,'fontweight','bold');
    if tt ==4
        [C2,h2]= contour(LON, LAT, SST_LS_trend_spec, [298.5:0.25:300.5],'w','linestyle','--','linewidth',1.1);
        clabel(C2,h2,[299:0.5:300],'color','w', 'labelspacing',400);
    end
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
    xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_' CloudType '_v1.jpg'],[0 0 10 8]);
    %xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_'
    %CloudType '_v2.jpg'],[0 0 10 8]);  with SST:

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



