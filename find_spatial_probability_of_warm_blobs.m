% This script is used to aggregate warm blob statistics sampled by RHB on different days
% and also blob statistics on different large scale condition dubbed
% (sugar-, gravel-, flowers-, fish-favored time);
% The goal is to compute: the innate likelihood of occurrence of
% warm anomalies at give location

% date: Jan 19, 2022; Feb 1, 2022 (updated)
% Note: 
clear all; clc; close all

%% load in all data:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
blobdatasvdir = [datadir filesep 'blob_data'];
product = 'g16';                  % 'gpblend'
TempAveTag = '5d-movingMean';     % '3d-movingMean'
SpatailFilterTag ='LPF600km'; % 'LPF600km';    % 'no_spatial_filtering'
cloudmaskTag = '';
if isempty(cloudmaskTag)
    casedir = strjoin({product, TempAveTag, SpatailFilterTag},'_');
else
    casedir = strjoin({product, TempAveTag, SpatailFilterTag, cloudmaskTag},'_');
end

figsvdir = [blobdatasvdir filesep casedir filesep 'figs'];


% datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
% datasvdir = [datadir filesep 'blob_data/detrended_large_scale_SSTgrad_updated'];

blob_files = dir([blobdatasvdir filesep casedir filesep  '*_individual_warm_blobs_mapSelCCthres0.5_blobSelCCthres0.1.mat']);
filenames = {blob_files.name};
nfiles = length(filenames);

%% load in all individual blob data:  (get the total number of sample times.
c0 = 1;
for i = 1:nfiles
    FN = filenames{i};
    load([blobdatasvdir filesep casedir filesep FN],'indiv_blobs_collection');
    
    nb = length(indiv_blobs_collection);
    cN = c0+nb-1;
    %try
    % temporarily skip the updated data.
    all_blobs(c0:cN) = indiv_blobs_collection;
    
    %all_blobstats(c0:cN) = indiv_blobs_stats;
    c0 = cN+1;
    % end
end

Nsamples = cN;      % total number of SST maps being sampled.

%% load cloudiness frequency estimate:
load([datadir filesep 'cloudiness_freq_ATOMIC.mat']);         % there is a local grid here: LON, LAT associated with the cloud frequency data.
P_clearsky_ng = 1-cloudy_freq_all;    % on the native grid;
% interpolation onto a coarser grid later.


%% for all individual blobs:
% 1. estimate probability of satellite being able to sense the signal of
% warm water:
BlobPixelLocs_all = get_pixelLocs_in_boundingbox(all_blobs);
%blob_centroids = vertcat(blobs_tmp.GeoLocs);
X = BlobPixelLocs_all(:,1);
Y = BlobPixelLocs_all(:,2);

locbin_size = 0.2;      % units: dgr
XBinEdge = [-59:0.2:-46];
YBinEdge = [8:0.2:20];

XBinCen = 0.5*(XBinEdge(1:end-1) + XBinEdge(2:end));
YBinCen = 0.5*(YBinEdge(1:end-1) + YBinEdge(2:end));

% pdf_all: over the 2 month period, the satellite see more warm anomalies
% in the yellow region. Doesn't imply the true location of the warm
% anomlies.
% just get the counts of warm anomlies events:
hc = histogram2(X,Y, XBinEdge, YBinEdge,'Normalization','count');
counts_all = hc.Values;

xres_dgr = 0.02;
P_detect_warmblobs = counts_all'./(Nsamples*(locbin_size/xres_dgr)^2);                 % on each map, the same pixel should only be counted once (a.k.a.: no overlap between individual blobs)

% find averaged cloud frequency on a the same grid box as the one showing the counts of warm anomalies
P_clearsky_smoothed = zeros(size(counts_all'));
[NY, NX] = size(LAT);
inc  = floor(locbin_size/xres_dgr);
for j = 1:size(counts_all,2)
    %yinds = (1+inc)*(j-1)+1:(1+inc)*j;
    condy = LAT>=YBinEdge(j)& LAT<YBinEdge(j+1);
    for i = 1:size(counts_all,1)
        %xinds = (1+inc)*(i-1)+1:(1+inc)*i;
        condx = LON>=XBinEdge(i)& LON<XBinEdge(i+1);
        mask = condx&condy;
        %         figure(9); clf;
        %         plot(LON,LAT,'.k');
        %         hold on
        %         plot(LON(mask), LAT(mask),'.r');
        %         hold off
        
        chunk = P_clearsky_ng(mask);
        
        P_clearsky_smoothed(j,i) = mean(chunk(:));
    end
end

if strcmp(cloudmaskTag, 'no_cloudmask')
    P_warmblobs_innate.all = P_detect_warmblobs;
else
    P_warmblobs_innate.all = P_detect_warmblobs./P_clearsky_smoothed;
end

figure(1); clf;
pcolor(XBinCen, YBinCen, P_warmblobs_innate.all);shading flat;
hb=colorbar;
set(get(hb,'xlabel'),'string','probability');
title('innate likelihood of existence of warm anormalies in Jan and Feb 2020')
hold on;
contour(XBinCen, YBinCen, P_warmblobs_innate.all,[0.1:0.1:0.6],'k');
caxis([0 0.8]);
xlabel('Longitude (^{\circ}E)');
ylabel('Latitude (^{\circ}N)');
xlim([-58, -48]);
ylim([8, 18]);  % or ylim([10, 20]
set(gca,'fontsize',14);
%xc_savefig(gcf,'Figs/likelihood_maps', ['innate_likelihood_of_warm_blobs_in_JanFeb2020_' casedir '.jpg'],[0 0 10 8]);
xc_savefig(gcf,figsvdir, ['innate_likelihood_of_warm_blobs_in_JanFeb2020_' casedir '_corrected.jpg'],[0 0 10 8]);

% the results are good.
% overlay increased cloud frequency contours of gravel on day time on to this likelihood map;
% overlay likelihood map onto the CF map will be better.
figure(2); clf;
pcolor(XBinCen, YBinCen, P_warmblobs_innate.all);shading flat; 
hb=colorbar;
set(get(hb,'xlabel'),'string','probability');
title('innate likelihood of existence of warm anormalies in Jan and Feb 2020')
hold on;
contour(XBinCen, YBinCen, P_warmblobs_innate.all,[0.1:0.1:0.6],'k');
% add the contours of cloudiness frequency:
[Ccld, hcld] = contour(LON, LAT, day.gravel.CF_anom, [0.2:0.2:0.6],'r');
clabel(Ccld, hcld, [0.2:0.2:0.6],'color','r');
caxis([0 0.8]);
xlabel('Longitude (^{\circ}E)');
ylabel('Latitude (^{\circ}N)');
xlim([-58, -48]);
ylim([8, 18]);  % or ylim([10, 20] 
set(gca,'fontsize',14);



%% group sampled blobs by U10 and LTS conditions:
% use a function to provide sorting.
grouped_blobs = group_blobs_by_U10_and_LTS_conditions(all_blobs);

% Then I can creat a map of 2D histogram that provide the count or
% probability of having a cloud blobs in the area bin. probability map.
cloud_types = {'gravel','flowers','sugar','fish'};

% should I actually find the probability for each type?
% The probability will be higher on the clear region obviously.
%% generating the spatial distribution of the anomalously warm pixels:
% check if the innate liklihood map of warm blob existence changes under
% different large scale atmospheric condition.
for tt = 1:4
    CloudType = cloud_types{tt};
    
    blobs_tmp = grouped_blobs.(CloudType);
    Nsample_type= length(blobs_tmp);
    % select day or night time blobs:
    blobs_time = [blobs_tmp.time];
   % thrs = (blobs_time - floor(blobs_time))*24;
%     nightTimeMask = false(size(blobs_time));
%     nightTimeMask(thrs<=10)=true;
%     nightTimeMask(thrs>21) =true;
%     idx_sel = find(nightTimeMask==1);
    
    blobs_stat_tmp = [blobs_tmp.stats_selected];
    %blob_EqvDiam = [blobs_stat_tmp.EqvDiam];
    
    % for each cloud type: do the following:
    BlobPixelLocs = get_pixelLocs_in_boundingbox(blobs_tmp);
    %blob_centroids = vertcat(blobs_tmp.GeoLocs);
    X = BlobPixelLocs(:,1);
    Y = BlobPixelLocs(:,2);
    
    figure(1);
    subplot(2,2,tt);
    % this probability haven't included information on size.
    hh = histogram2(X,Y, XBinEdge, YBinEdge,'Normalization','count',...
        'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
    colorbar
    title(CloudType);
    %caxis([0 0.03]);
    
    %% estimate the probability of detecting warm anomalies:
    counts_type=hh.Values';
    P_detect_warmblobs_type = counts_type./(Nsample_type*(locbin_size/xres_dgr)^2);                 % on each map, the same pixel should only be counted once (a.k.a.: no overlap between individual blobs)
    
    P_clearsky_type = 1-cloudy_freq.(CloudType);
    P_clearsky_smoothed_type = zeros(size(counts_type));
    [NY, NX] = size(LAT);
    inc  = floor(locbin_size/xres_dgr);
    for j = 1:size(counts_all,2)
        %yinds = (1+inc)*(j-1)+1:(1+inc)*j;
        condy = LAT>=YBinEdge(j)& LAT<YBinEdge(j+1);
        for i = 1:size(counts_all,1)
            %xinds = (1+inc)*(i-1)+1:(1+inc)*i;
            condx = LON>=XBinEdge(i)& LON<XBinEdge(i+1);
            mask = condx&condy;
            %         figure(9); clf;
            %         plot(LON,LAT,'.k');
            %         hold on
            %         plot(LON(mask), LAT(mask),'.r');
            %         hold off
            
            chunk = P_clearsky_type(mask);
            
            P_clearsky_smoothed_type(j,i) = mean(chunk(:));
        end
    end
    
    if strcmp(cloudmaskTag, 'no_cloudmask')
        P_warmblobs_innate.(CloudType) = P_detect_warmblobs_type;
    else
        P_warmblobs_innate.(CloudType) = P_detect_warmblobs_type./P_clearsky_smoothed_type;
    end
    
    hfig2= figure(2); %clf;
    subplot(2,2,tt);
    pcolor(XBinCen, YBinCen, P_warmblobs_innate.(CloudType));shading flat;
    hb=colorbar;
    set(get(hb,'xlabel'),'string','probability');
    title({'innate likelihood of existence of warm anormalies'; [CloudType '-favored conditions']})
    hold on;
    [C,hl]=contour(XBinCen, YBinCen, P_warmblobs_innate.(CloudType), [0.1:0.1:0.6],'k');
    clabel(C,hl, [0.2:0.2:0.8],'color','k');
    caxis([0 0.8]);
    xlabel('Longitude (^{\circ}E)');
    ylabel('Latitude (^{\circ}N)');
    xlim([-58, -48]);
    ylim([8, 18]);  % or ylim([10, 20]
    set(gca,'fontsize',14);
    
    
    
end
xc_savefig(hfig2, figsvdir, ['likelihood_of_warmblob_existence_4types_' casedir '.jpg'],[0 0 12 10]);

save([blobdatasvdir filesep casedir filesep 'likelihood_of_warmblob.mat'],'XBinCen','YBinCen','P_warmblobs_innate');



%% look at warm blob statistics in the region of interest:
gravel_blobs = grouped_blobs.gravel;
ROI(1).coverage = [-55, -51; 13.5, 15.5];
ROI(1).name = 'RHB_sampled_Hotspots';

ROI(2).coverage = [-56, -54; 10, 12];
ROI(2).name = 'non_hotspot_area';

for ig = 1 %:2
    examine_region = ROI(ig).coverage;
    region_name = ROI(ig).name;
    
    % examine_region = [-56, -54; 10, 12];
    % region_name = 'non_hotspot_area';
    
    % get locations of all the blobs;
    BlobCenLocs = vertcat(gravel_blobs.GeoLocs);
    for i = 1:length(gravel_blobs)
        gravel_blobs(i).timevec = gravel_blobs(i).time.*ones(size(gravel_blobs(i).max_SSTa));
    end
    BlobTime = [gravel_blobs.timevec];
    bx = BlobCenLocs(:,1); by = BlobCenLocs(:,2);
    
    xres = 2;  % 2km
    gravel_blob_stats = get_blob_characters(gravel_blobs, xres);
    
    % select blobs within the examination region:
    rvx = [examine_region(1,1), examine_region(1,1), examine_region(1,2), examine_region(1,2), examine_region(1,1)];
    rvy = [examine_region(2,1),examine_region(2,2), examine_region(2,2), examine_region(2,1), examine_region(2,1)];
    [inflag, onflag] = inpolygon(bx,by, rvx, rvy);
    selflag = inflag|onflag;
    
    varlist = fieldnames(gravel_blob_stats);
    for iv = 1:length(varlist)
        VN = varlist{iv};
        selected_blobs.(VN) = gravel_blob_stats.(VN)(selflag);
    end
    selected_blobs.time = BlobTime(selflag);
    
    % separte them into blobs detected in daytime and detected in night time;
    thrs = (selected_blobs.time - floor(selected_blobs.time))*24;
    nightTimeMask = false(size(thrs));
    nightTimeMask(thrs<=10)=true;    % locally: 6pm->6am next day.
    nightTimeMask(thrs>21) =true;
    
    idx_night = find(nightTimeMask==1);
    idx_day = find(nightTimeMask==0);
    
    %% plot the distribution of statistics or something else?
    % 1. I wanted to know if there are more night time blobs than day time
    % blobs (nope.); about the same, 51.4% sampled during the day and 48.6% catched at
    % night.
    
    % look at size, SST anomaly and background SST, AxisRatio, Orientation for
    % day time and night time blobs;
    
    plotvars = {'SST_anom','EqvDiam_km','AxisRatio','SST_bg','Orientation'};
    xlabels ={'max. SST anomaly (^{\circ}C)', 'Equivalent Diameter (km)', 'Axis Ratio', ...
        'ave. background SST (^{\circ}C)', 'Orientation'};
    edges.SST_anom = [0:0.1:1];
    edges.EqvDiam_km = [0:20:200];     % km
    edges.AxisRatio = [1:0.2:4];
    edges.SST_bg = [298:0.2:301]-273;
    edges.Orientation = [-90:10:90];
    RHBcolor = [0.8500 0.3250 0.0980];
    
    xticks.SST_anom = [0:0.2:1];
    xticks.EqvDiam_km = [0:40:200];
    xticks.AxisRatio = [1:0.5:4];
    xticks.SST_bg = [25:0.5:28];
    xticks.Orientation = [-90:30:90];
    
    yrange.SST_anom = [0 250];
    yrange.EqvDiam_km=[0 350];
    yrange.AxisRatio = [0 120];
    yrange.SST_bg = [0 300];
    yrange.Orientation = [0 120];
    
    % distribution plot:
    figure(1); clf; clear h
    for iv = 1:length(plotvars)
        VN = plotvars{iv};
        val_night = selected_blobs.(VN)(idx_night);
        val_day = selected_blobs.(VN)(idx_day);
        
        if iv==4
            val_night = val_night-273.15;
            val_day = val_day-273.15;
        end
        subplot(2,3,iv);
        h(1) = histogram(val_night, edges.(VN));
        hold on
        h(2) =histogram(val_day, edges.(VN)); %, 'FaceAlpha', 0.4);
        %yrange = get(gca,'ylim');
        %     plot(median(val_night)*ones(1,2),yrange,'color', 'c','linewidth',1.2);
        %     plot(median(val_day)*ones(1,2),yrange,'color', 'm','linewidth',1.2);
        
        xlabel(xlabels{iv});
        ylabel('count');
        %title(VN);
        if iv == 1
            legend(h, {'nighttime','daytime'});
            
        end
        set(gca,'fontsize',13);
        set(gca,'xtick', xticks.(VN),'TickDir','both');
        %ylim(yrange.(VN));
        grid on
        
    end
    subplot(2,3,6)
    labels = {['night-time (' num2str(length(idx_night)) ')'], ...
        ['daytime (' num2str(length(idx_day)) ')']};
    tp  = pie([length(idx_night), length(idx_day)]);
    tp(2).FontSize = 13;
    %tp(2).BackgroundColor = h(1).Cont
    tp(4).FontSize= 13;
    lgd = legend(labels);
    lgd.Location = 'southoutside';
    xc_savefig(gcf,figsvdir,['gravel_' region_name '_blobstats_dayVSnight_autoyrange.jpg'],[0 0 12 8]);
    
    
    %% plot scatter plots using only the size, strength and color-coded by background SST;  % remove orientation and axis-ratio.
    % other plot, statistical distribution:
    % any differences betweent the day time and night time blobs?
    figure(3); clf;
    for i = 1:2
        subplot(2,1,i)
        if i==1
            idx_sel = idx_night;
        else
            idx_sel = idx_day;
        end
        
        SSTa = [selected_blobs.SST_anom(idx_sel)];
        AxisRatio = [selected_blobs.AxisRatio(idx_sel)];
        EqvDiam = [selected_blobs.EqvDiam_km(idx_sel)];
        SSTbg = [selected_blobs.SST_bg(idx_sel)]-273.15;
        Orientation = [selected_blobs.Orientation(idx_sel)];
        
        %% plot data where the EqvDiam is less than 75th percentile.
        sizemask = EqvDiam<prctile(EqvDiam, 75);
        %scatter(SSTa(sizemask),  EqvDiam(sizemask), 40, SSTbg(sizemask), 'filled','MarkerEdgeColor','k', 'MarkerFaceAlpha',1);
        scatter(SSTa,  EqvDiam, EqvDiam*1.1, SSTbg, 'filled','MarkerEdgeColor','k', 'MarkerFaceAlpha',1);

        hold on
        %% plot data where EqvDiam is larger than 75th percentile:
        %scatter(Orientation(~sizemask), SSTa(~sizemask),  EqvDiam(~sizemask)*1.1, SSTbg(~sizemask), 'filled','MarkerEdgeColor','k', 'MarkerFaceAlpha',1.0);
        % this is actually easier to see the differences.
        

        stats =find_statistical_measures([SSTa; EqvDiam]);
        plot(stats(1).median, stats(2).median, '*m', 'markersize',20,'linewidth',2);
        plot(stats(1).mean, stats(2).mean,'dm', 'markersize',15,'linewidth',2);
        hold on;
        IQRbox.x = [stats(1).InterQuartile(1), stats(1).InterQuartile(1), stats(1).InterQuartile(2), stats(1).InterQuartile(2), stats(1).InterQuartile(1)];
        IQRbox.y =[stats(2).InterQuartile(1), stats(2).InterQuartile(2), stats(2).InterQuartile(2), stats(2).InterQuartile(1),stats(2).InterQuartile(1)];
        plot(IQRbox.x, IQRbox.y, '-k','linewidth',1.2);
        
        leftHand.x = [stats(1).lowerWhisker, stats(1).InterQuartile(1)];
        leftHand.y = [stats(2).median, stats(2).median];
        rightHand.x = [stats(1).upperWhisker, stats(1).InterQuartile(2)];
        rightHand.y = leftHand.y;
        leftEdge.x = [stats(1).lowerWhisker,stats(1).lowerWhisker];
        rightEdge.x = [stats(1).upperWhisker,stats(1).upperWhisker];
        vEdge.y = [stats(2).median*1.1, stats(2).median*0.9];
        
        plot(leftHand.x, leftHand.y,'--k','linewidth',1.2); plot(rightHand.x, rightHand.y,'--k','linewidth',1.2);
        plot(leftEdge.x, vEdge.y,'-k','linewidth',1.2); plot(rightEdge.x, vEdge.y,'-k','linewidth',1.2);
        
        % create upper and lower edges:
        lowerWhisk.y = [stats(2).lowerWhisker, stats(2).InterQuartile(1)];
        Whisk.x = [stats(1).median, stats(1).median];
        upperWhisk.y = [stats(2).upperWhisker, stats(2).InterQuartile(2)];
        hEdge.x = [stats(1).median*1.1, stats(1).median*0.9];
        lowerEdge.y = [stats(2).lowerWhisker, stats(2).lowerWhisker];
        upperEdge.y = [stats(2).upperWhisker, stats(2).upperWhisker];
        
        plot(Whisk.x, lowerWhisk.y, '--k','linewidth',1.2); plot(Whisk.x, upperWhisk.y, '--k','linewidth',1.2);
        plot(hEdge.x, lowerEdge.y, '--k','linewidth',1.2); plot(hEdge.x, upperEdge.y, '--k','linewidth',1.2);
        
        % plot threshold:
        plot([0.37, 0.37],[0, 180],'--r','linewidth',1.2);
        hold on;
        plot([0, 0.8],[38, 38], '--r','linewidth',1.2);

        
        % tag blobs with dates:
        %timetag = cellstr(datestr(blob_time,'mmmdd'));
        %text(SSTa, AxisRatio, timetag,'fontsize',11);
        hb = colorbar;
        colormap(turbo);
        set(get(hb,'xlabel'),'string','background SST (^{\circ}C)');
        xlabel('\DeltaSST_{max} (^{\circ}C)');
        ylabel('Equivalent Diameter (km)');
        set(gca,'fontsize',16, 'xtick',[0:0.2:0.8]);
        grid on
        ylim([0 180]);
        xlim([0, 0.8]);
        caxis([25.5, 27.5]);
        title(['Region ' char('A'+ig-1) ' in ' labels{i}]);
        
                % add a subaxes to show the scale of 75th EqvDiam:
%         ax1 = gca;
%         ax1_pos = ax1.Position;
%         ax2= axes('Position', [ax1_pos(1)+0.8*ax1_pos(3), ax1_pos(2)+0.85*ax1_pos(4), 0.15*ax1_pos(3), 0.15*ax1_pos(4)]);
%         
%         hs(1) = scatter(ax2, 0.5,2,prctile(EqvDiam, 75)*1.1,'k','linewidth',1.2,'MarkerFaceColor',[0.5 0.5 0.5]);
%         text(0.5,1, ['75th EqvDiam = ' num2str(round(prctile(EqvDiam, 75))) 'km'],'HorizontalAlignment','center');
%         ax2.XAxis.Color = 'none';
%         ax2.YAxis.Color='none';
        %ylim([0 1.5])
    end
    % perhaps plot the median value of both SSTa and Orientation and its 25th,
    % 75th percentile; as well as standard deviation of both number.
    %varnames={'SSTa','Orientation'};
    xc_savefig(gcf,figsvdir,[region_name '_dayVSnight_BlobStatsScatters_v2.jpg'],[0 0 10 10]);
end

%% obsolete:
if 1==0
    figure(6);
    % this probability haven't included information on size.
    h = histogram2(X,Y, XBinEdge, YBinEdge,'Normalization','pdf', ...
        'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
    colorbar
    title('Probability of warm blob detection from the 2-month data');
    
    pdf_all = h.Values;
    
    hold on;
    [C,h]=contour(XBinCen, YBinCen, pdf_blobdetection',[0.005:0.0025:0.03],'w');
    
    save([datadir filesep 'probability_density_function_for_warmblob_detection_ATOMIC.mat'],...
        'pdf_conditioned','pdf_all','XBinCen','YBinCen','XBinEdge','YBinEdge');
    
    %% Jan 25, 2022: get warm blob characteristics for different large scale conditions:
    %  - Note: for each conditions, I can further separate the statistics for day
    %  time and night time (to remove potential influence from diurnal warming)
    %
    cloud_types = {'gravel','flowers','sugar','fish'};
    res = 2;    % units: km
    % Note: add averaged cloud percentages for each type (cloud contamination
    % percentage in each SST map grouped into the type)
    figure(3); clf;
    for tt = 1:4
        CloudType = cloud_types{tt};
        
        blobs_tmp = grouped_blobs.(CloudType);
        stats = get_blob_characters(blobs_tmp, res);
        
        % for each cloud type: do the following:
        subplot(2,2,tt)
        % this probability haven't included information on size.
        % not so informative, need to plot it in a different way?
        %% scatter plots:
        scatter(stats.SST_bg-273.15, stats.EqvDiam_km, 30, stats.SST_anom, 'filled', 'marker','o');
        hold on;
        plot(mean(stats.SST_bg)-273.15*ones(1,2), [0 700],'--r','linewidth',1.2);
        xlabel('background SST for the warm anoamlies');
        ylabel('Eqivalent Diameter (km)');
        %plot(median(stats.SST_bg)-273.15*ones(1,2), [0 700],'--m','linewidth',1.2);
        %plot(median(stats.SST_anom), median(stats.EqvDiam_km),'+w','linewidth',3,'markersize',10);
        %% histogram
        %     histogram2(stats.SST_anom, stats.EqvDiam_km, [0:0.1:2.5], [0:20:600],'Normalization','pdf', ...
        %         'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
        colorbar;
        caxis([0, 1]);
        %caxis([26,28]);
        xlim([24.5, 28.5]);
        ylim([0 700]);
        title(CloudType);
        
    end
    
    % I think that I can check where the anomalies warmer than certin
    % threshold locate in space. as well as anomalies with large size.  the spatial
    % distribution map can be generate for warm blobs that satisfy a particular
    % criteria.
    select_warm_blobs_with_criteria();    % time (day / night)
    % SST anomaly thresholds
    % size
    location_distribution_of_selected_warm_blobs();
    
    end