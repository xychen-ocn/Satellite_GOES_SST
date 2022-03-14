clear all; clc; close all;

dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
datasvdir = [dataroot filesep 'blob_data/L4test'];
figsvdir = [datasvdir filesep 'figs'];

matFN = 'L4_gpblend_daily_WarmSSTanom_blobs_and_cloudiness_larger_feature_domain.mat';
load([datasvdir filesep matFN]);

% composite blobs by large scale atmospheric conditions:
for i = 1:length(indiv_blobs_collection)
    tmp_struct = indiv_blobs_collection(i);
    tmp_struct.cloudInfo = cloudiness_downwind_FCN{i};
    all_blobs(i) = tmp_struct;
    clear tmp_struct
end

grouped_blobs = group_blobs_by_U10_and_LTS_conditions(all_blobs);    % this may need to be improved to use daily mean large scale conditions.

tmp = [all_blobs.cloudInfo];
all_windspd =[tmp.mean_wspd];
histogram(all_windspd, 20);

% I can further constrain locations of the warm anomalies.
search_region = [-59, -48; 8, 18];
svx = [-59, -59, -48, -48, -59];
svy = [8, 18, 18, 8, 8];

xlabelstr = {'crxwind'}; %X_{normalized}
ylabelstr = {'downwind'}; %Y_{normalized}


%% all cloudiness:
allblob_locs= vertcat(all_blobs.GeoLocs);

inmask = inpolygon(allblob_locs(:,1), allblob_locs(:,2), svx, svy);
selids = find(inmask==1);

all_cloudInfo = [all_blobs.cloudInfo];

all_CFdata = grid_data_into_windaligned_coord(all_cloudInfo(selids));
CF_allcomp_ave = mean(all_CFdata.cloudyfreq,3,'omitnan');
%CF_allcomp_cc = sum(all_CFdata.cloudycnt,3,'omitnan')./sum(all_CFdata.samplesz);

figure(1);
% mean value:
subplot(1,2,1);
pcolor(all_CFdata.XX, all_CFdata.YY, CF_allcomp_ave);shading flat;
hold on
circle(0,0, 1,'k',1.3);
hb=colorbar; set(get(hb,'xlabel'), 'string','cloudiness');
colormap(hsub, turbo)
caxis([0.1, 0.6])
hold on
[hl_all, c_all]=contour(all_CFdata.XX, all_CFdata.YY, CF_allcomp_ave,[0.3:0.01:0.6],'k');
clabel(hl_all, c_all, [0.3:0.01:0.6],'color','k');
set(gca,'fontsize',12);
title(['cloudiness ALL']);
xlabel(xlabelstr);
ylabel(ylabelstr);
axis('square');

% remove domain mean value:
hsub=subplot(1,2,2);
domain_mean = mean(CF_allcomp_ave,'all');
pcolor(all_CFdata.XX, all_CFdata.YY, CF_allcomp_ave-domain_mean);shading flat;
hold on
circle(0,0, 1,'k',1.3);
hb=colorbar; set(get(hb,'xlabel'), 'string','anomaly');
colormap(hsub, redblue)
caxis([-0.05, 0.05])
hold on
[hl_all, c_all]=contour(all_CFdata.XX, all_CFdata.YY, CF_allcomp_ave-domain_mean,[-0.05:0.01:0.05],'k');
clabel(hl_all, c_all, [-0.05:0.01:0.05],'color','k');
set(gca,'fontsize',12);
xlabel(xlabelstr);
%ylabel(ylabelstr);
title(['- domain mean (' num2str(domain_mean,'%5.3f') ')']);
axis('square');




    
%% cloudiness in different atmospheric regimes:
cloudtypes={'flowers','gravel','sugar','fish'};
gridX = [-5:0.1:5]; gridY = gridX;

for i = 1:4
    CT = cloudtypes{i};
    allblob_locs = vertcat(grouped_blobs.(CT).GeoLocs);
    
    inmask = inpolygon(allblob_locs(:,1), allblob_locs(:,2), svx, svy);
    selids = find(inmask==1);
    
    all_cloudInfo = [grouped_blobs.(CT).cloudInfo];
    
    
   %% compute cloudiness in daytime and in nighttime separately. 
   for j = 1:length(selids)
       k = selids(j);
       SST_cmasks = all_cloudInfo(k).SST_cloudmask;
       [cldmask_day_tmp, cldmask_night_tmp] = group_data_into_daytime_and_nighttime_samples(SST_cmasks, all_cloudInfo(k).time);
       % compute cloudiness separately:
       cldmask_day(j).cloudfreq=compute_cloudfreq(cldmask_day_tmp.SSTmaps);
       cldmask_night(j).cloudfreq=compute_cloudfreq(cldmask_night_tmp.SSTmaps);
       cldmask_day(j).cldmasks = cldmask_day_tmp.SSTmaps;
       cldmask_night(j).cldmasks = cldmask_night_tmp.SSTmaps;
       
       cldmask_day(j).WindAligned_Coord = all_cloudInfo(k).WindAligned_Coord;
       cldmask_night(j).WindAligned_Coord = all_cloudInfo(k).WindAligned_Coord;
       
   end
    
    all_CFdata = grid_data_into_windaligned_coord(all_cloudInfo(selids), gridX, gridY);
    CF_RegimeComp_ave.(CT) = mean(all_CFdata.cloudyfreq,3,'omitnan');
    %CF_allcomp_cc = sum(all_CFdata.cloudycnt,3,'omitnan')./sum(all_CFdata.samplesz);
    
    daytime_CFdata = grid_data_into_windaligned_coord(cldmask_day, gridX, gridY);
    CF_RegimeDaytimeComp_ave.(CT) = mean(daytime_CFdata.cloudyfreq,3,'omitnan');
    
    nighttime_CFdata = grid_data_into_windaligned_coord(cldmask_night, gridX, gridY);
    CF_RegimeNgttimeComp_ave.(CT) = mean(nighttime_CFdata.cloudyfreq,3,'omitnan');


    
end
    
figure;
pos =customize_subplot_size(4, 3, 0.15,0.04);
dataIn = CF_RegimeComp_ave;
for i  = 1 :4
    CT = cloudtypes{i};
    
    xlabelstr = {'crxwind'}; %X_{normalized}
    ylabelstr = {'downwind'}; %Y_{normalized}
    hsub = subplot(4,3,1+3*(i-1));
    pcolor(all_CFdata.XX, all_CFdata.YY, dataIn.(CT));shading flat;
    hold on
    circle(0,0, 1,'k',2);
    hb=colorbar; set(get(hb,'xlabel'), 'string','cloudiness');
    colormap(hsub, turbo)
    caxis([0.1, 0.6])
    hold on
    [hl_all, c_all]=contour(all_CFdata.XX, all_CFdata.YY, dataIn.(CT),[0.3:0.01:0.6],'k');
    clabel(hl_all, c_all, [0.3:0.01:0.6],'color','k');
    set(gca,'fontsize',12);
    title(['cloudiness: ' CT]);
    if i ==4
        xlabel(xlabelstr);
    end
    ylabel(ylabelstr);
    axis('square');
    set(hsub, 'pos',pos{1+3*(i-1)});
    
    hsub = subplot(4,3,2+3*(i-1));
    domain_mean = mean(dataIn.(CT),'all');
    pcolor(all_CFdata.XX, all_CFdata.YY, dataIn.(CT)-domain_mean);shading flat;
    hold on
    circle(0,0, 1,'k',2);
    hb=colorbar; set(get(hb,'xlabel'), 'string','anomaly');
    colormap(hsub, redblue)
    caxis([-0.1, 0.1])
    hold on
    [hl_all, c_all]=contour(all_CFdata.XX, all_CFdata.YY, dataIn.(CT)-domain_mean,[-0.08:0.02:0.08],'k');
    clabel(hl_all, c_all, [-0.08:0.02:0.08],'color','k');
    set(gca,'fontsize',12);
    if i==4
        xlabel(xlabelstr);
    end
    %ylabel(ylabelstr);
    title(['- domain mean (' num2str(domain_mean,'%5.3f') ')']);
    axis('square');
    set(hsub, 'pos',pos{2+3*(i-1)});
    
    hsub = subplot(4,3,3+3*(i-1));
    downwind_mean = mean(dataIn.(CT));
    pcolor(all_CFdata.XX, all_CFdata.YY, dataIn.(CT)-downwind_mean);shading flat;
    hold on
    circle(0,0, 1,'k',2);
    hb=colorbar; set(get(hb,'xlabel'), 'string','anomaly');
    colormap(hsub, redblue)
    caxis([-0.1, 0.1])
    hold on
    [hl_all, c_all]=contour(all_CFdata.XX, all_CFdata.YY, dataIn.(CT)-downwind_mean,[-0.08:0.02:0.08],'k');
    clabel(hl_all, c_all, [-0.08:0.02:0.08],'color','k');
    set(gca,'fontsize',12);
    if i==4
        xlabel(xlabelstr);
    end
    %ylabel(ylabelstr);
    title('- downwind-mean');
    axis('square');
    set(hsub, 'pos',pos{3+3*(i-1)});
    
end
xc_savefig(gcf, figsvdir, 'averaged_daily_cloudiness_4CloudTypes_5x5domain.jpg',[0 0 12 10]);




%% cloudiness composite based on wind speed and background SST:
% do composite of cloudiness based on wind speed, background SST;
% try gravel first:
for i  = 1 %:4
    CT = cloudtypes{i};
    
   % CT = 'gravel';
    bloblocs = vertcat(grouped_blobs.(CT).GeoLocs);
    
    inmask = inpolygon(bloblocs(:,1), bloblocs(:,2), svx, svy);
    selids = find(inmask==1);
    
    tmp = [grouped_blobs.(CT).cloudInfo];
    gravel_windspd = [tmp.mean_wspd];
    gravel_SSTbg = [grouped_blobs.(CT).ave_SSTbg]-273.15;
    gravel_SSTmax = [grouped_blobs.(CT).max_SSTa];
    blob_stats = [grouped_blobs.(CT).stats_selected];
    gravel_blobsize = [blob_stats.EqvDiam]*5.5;  % km
    
    figure(11);
    subplot(4,1,1)
    histogram(gravel_windspd(selids));%,[6:1:14]);
    title('wind speed');
    
    subplot(4,1,2)
    histogram(gravel_SSTbg(selids));%,[24:1:28]);
    title('large scale background SST');
    
    subplot(4,1,3)
    histogram(gravel_SSTmax(selids));%,[0:0.1:0.6]);
    title('max \DeltaSST');
    
    subplot(4,1,4)
    histogram(gravel_blobsize(selids)); %,[0:20:300]);
    title('feature EqvDiam (km)');
    
    %% check distribution (2D histogram) of SSTbg and max delta SST:
    figure(10);
    subplot(3,1,1)
    h = histogram2(gravel_SSTbg(selids),gravel_SSTmax(selids), [25:0.5:28],[0.1:0.1:0.5],...
        'DisplayStyle','tile','ShowEmptyBins','on');
    hb=colorbar; set(get(hb,'xlabel'),'String','occurence');
    set(gca,'ytick',[0.1:0.1:0.5]);
    
    xlabel('background large-scale SST (Celsius)');
    ylabel('\DeltaSST (Celcius)');
    set(gca,'fontsize',14);
    
    subplot(3,1,2)
    h = histogram2(gravel_SSTbg(selids),gravel_blobsize(selids), [25:0.5:28],[0:20:200],...
        'DisplayStyle','tile','ShowEmptyBins','on');
    hb=colorbar; set(get(hb,'xlabel'),'String','occurence');
    set(gca,'ytick',[0:20:200]);
    xlabel('background large-scale SST (Celsius)');
    ylabel('EqvDiam (km)');
    set(gca,'fontsize',14);
    
    subplot(3,1,3)
    h = histogram2(gravel_SSTmax(selids),gravel_blobsize(selids), [0.1:0.1:0.5],[0:20:200],...
        'DisplayStyle','tile','ShowEmptyBins','on');
    hb=colorbar; set(get(hb,'xlabel'),'String','occurence');
    set(gca,'ytick',[0:20:200]);
    xlabel('\DeltaSST (Celsius)');
    ylabel('EqvDiam (km)');
    set(gca,'fontsize',14);

    
    
    % simply check the cloudiness composite for wind speed;
    % for each wind speed group, we can further check bg_SST: (to-do
    % tomorrow);
    if strcmp(CT,'gravel')
        wspd_edges = [6:2:14];  %gravel
        %SSTbg_edges = [25.5:0.5:27.5];
    elseif strcmp(CT,'fish')
        wspd_edges = [0:2:8];    %fish
    elseif strcmp(CT,'sugar')
        wspd_edges =[3:1:8];     %sugar
    end
    
    %
    SSTbg_edges = [25.5:0.5:28];
    
    SSTa_edges = [0.1:0.1:0.5];
    size_edges = [0:20:100];
    
    matrixname = 'SSTbg';
    matrix = gravel_SSTbg(selids);
    matrix_edges = SSTbg_edges;
    
%     matrixname = 'SSTa';
%     matrix = gravel_SSTmax(selids);
%     matrix_edges = SSTa_edges;
%     
    
    
    clear gids sstids
    Widx=discretize(gravel_windspd(selids), wspd_edges);
    for id = 1:length(wspd_edges)-1
        % loop through different bins of wind speed
        gids{id} = find(Widx==id);
        
        % for each wind speed bin, sort blobs further by background SST:
        %SSTidx = discretize(gravel_SSTmax(gids{id}), SSTa_edges);
        SSTidx = discretize(matrix(gids{id}), matrix_edges);
        
        for k = min(SSTidx):max(SSTidx)
            tmp = find(SSTidx==k);
            if ~isempty(tmp)
                sstids{id, k} = tmp;
            end
        end
        %
    end
    
    % use the index above to composite cloudiness data:
    % what to do?
    dataIn = [grouped_blobs.(CT).cloudInfo];
    [CF_composites, coord, samplesize, hfigs] = composite_data_in_windaligned_coord(dataIn(selids), sstids,'aveCF');
    figsvdir = [dataroot filesep 'blob_data/L4test/figs'];
    if ~exist(figsvdir)
        mkdir(figsvdir)
    end
    figname_prefix = [CT  '_ATOMIC_cloudiness_downwindCFN_sortby_windspeed4bins_and_' matrixname ];
    for i =1:length(hfigs)
        figname = [figname_prefix '_' num2str(i,'%2.2i') '_largerDomain.jpg'];
        xc_savefig(hfigs(i), figsvdir, figname, [0 0 12 10]);
    end
    
end



%% explore diunral variations:



% plot results:
 %% put results from all conditions;
 type = 'domain_mean';
 hfigs = visualize_composite(coord, CF_composites, type, samplesize);
 figname_prefix = [CT  '_ATOMIC_cloudiness_downwindCFN_sortby_windspeed4bins_and_' matrixname ];
for i =1:length(hfigs)
    figname = [figname_prefix '_rmved_domain_mean.jpg'];
    xc_savefig(hfigs(i), figsvdir, figname, [0 0 12 10]);
end


 type = 'downwind_mean';
 hfigs = visualize_composite(coord, CF_composites, type, samplesize);

 


% the cloudiness map: need to map on the gridded wind_aligned coord for
% composite.  (to-do tomorrow)

    