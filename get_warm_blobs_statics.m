% This script will summarize the basic statistics of warm blobs in the
% ATOMIC broder area ()
% this should be a function
%
% Note: the same "warm blob" can be detected at different time steps due to 
%       different reasons:
%     a. for given dT_thres, the warm blob presents for a certain amount of
%     time due to passing clouds, also, the warm blob can be moving around.
%     b. At given time, as dT_thres increases, the same warm blob can be
%     detected wil changing size due to a threshold change.
% 
%  If we assume that the warm blob doesn't move much from hours to hours (quasi-stationary), then we
%  identify the individual warm blob by its location in space. 
%  At given time, we will take the largest warm blob from the warm blob doppelgangers 
%  as the size and shape of the true warm blob.
%
%  Then, for each warm blob, we can get the daily-averaged maximum SST
%  anomalies (to account for its magnitude change in time.)
%
%  Then we can basically obtain the 3S: Size, Shape, and Strength.
%

% This script will be used to get the RHB sampled blobs' statistics;

%% load in data for a day of interest as an example to build up the scripts:
clear all; clc; close all;
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
datasvdir = [datadir filesep 'blob_data/detrended_large_scale_SSTgrad_updated'];

matfiles = dir([datasvdir filesep '*RHB_sampled_blobs.mat']);
matFNs = {matfiles.name};
t0 = 1;
for i = 1:length(matfiles)
    load([datasvdir filesep matFNs{i}]);   
    nb = length(RHB_sampled_blobs);
    tN = t0+nb-1;
    all_RHBblobs(t0:tN) = RHB_sampled_blobs;
    all_RHBblobs_stats(i) = RHB_blobs_stats;
    t0 = tN+1;
end

load([datadir filesep 'RHB_wavelet_coherent_segment_Info.mat']);
RHB_days = unique(cellstr(datestr(RHB_segtime,'mmmdd')));                                 % include days within the segment time.

for i = 1:length(matfiles)
    tmp = strsplit(matFNs{i},'_');
    
    RHB_catch_dateID{i} = tmp{1};
end

format = 'yyyy-mm-ddTHH:MM:SS';
notSWC = {['2020-02-03T19:00:00'; '2020-02-04T19:30:00']; ...
          ['2020-01-23T22:00:00'; '2020-01-24T18:05:00'];
          ['2020-02-12T15:00:00'; '2020-02-12T23:59:00']};
isSWC = {['2020-01-09T00:00:00'; '2020-01-10T01:20:00']; ...
          ['2020-01-22T20:00:00'; '2020-01-23T10:30:00']; ...
          ['2020-02-06T18:00:00'; '2020-02-07T23:00:00']; ...
          ['2020-02-12T06:00:00'; '2020-02-12T14:30:00']};

for i = 1:length(notSWC)
    notSWC_timerange(:,i) = datenum(notSWC{i},format);
end

for i = 1:length(isSWC)
    isSWC_timerange(:,i) = datenum(isSWC{i},format);
end

%%
% find stats within the notSWC and isSWC time range;

RHBblobs.isSWC = group_RHBblobs_by_WaveletCoherenceflag(isSWC_timerange, all_RHBblobs);
RHBblobs.notSWC = group_RHBblobs_by_WaveletCoherenceflag(notSWC_timerange, all_RHBblobs);

flagNames = fields(RHBblobs);

%% 2. get the 3S from these RHB sampled warm blobs. 
stats.Strength = [all_RHBblobs_stats.SST_anom];
stats.Size = [all_RHBblobs_stats.EqvDiam_km];
stats.Shape = [all_RHBblobs_stats.AxisRatio];
stats.backgroundSST = [all_RHBblobs_stats.SST_bg] - 273;
stats.MajorAxis = [all_RHBblobs_stats.MajAxisLen];
stats.MinorAxis = [all_RHBblobs_stats.MinAxisLen];

edges.Strength = [0:0.1:1];
edges.Size = [0:20:200];     % km
edges.Shape = [1:0.2:4];
edges.backgroundSST = [299:0.2:301]-273;
edges.MajorAxis = [0:20:200];
edges.MinorAxis = edges.MajorAxis;
RHBcolor = [0.8500 0.3250 0.0980];

varlist = fieldnames(stats);
xlabels ={'SST anomaly (K)', 'Equivalent Diameter (km)', 'Axis Ratio', ...
    'ave. background SST (^{\circ}C)', 'Major Axis Length (km)','Minor Axis Length (km)'};

figure(1);
for iv = 1:length(varlist)
    VN = varlist{iv};
    val = stats.(VN);
    subplot(2,3,iv);
    histogram(val, edges.(VN));
    hold on
    yrange = get(gca,'ylim');
    plot(median(val)*ones(1,2),yrange,'color', RHBcolor,'linewidth',1.2);

    xlabel(xlabels{iv});
    ylabel('count');
    title(VN);
end


%%
figure(2); clf;
marker= {'o','d'};
for i = 1:2
    flag = flagNames{i};
    SSTa = [RHBblobs.(flag).max_SSTa];
    tmp_stats =  [RHBblobs.(flag).stats_selected];
    AxisRatio = unique([tmp_stats.MajAxisLen])./[tmp_stats.MinAxisLen];
    EqvDiam = [tmp_stats.EqvDiam];
    SSTbg = [RHBblobs.(flag).ave_SSTbg]-273.15;
    blob_time = [RHBblobs.(flag).time];
    
    scatter(SSTa,AxisRatio, EqvDiam*1.1, SSTbg, 'filled','marker',marker{i});
    if i==2
        scatter([SSTa, NaN],[AxisRatio, NaN], [EqvDiam*1.1, NaN], [SSTbg, NaN], 'filled','marker',marker{i}, ...
            'markerEdgeColor','k');
    end
    hold on;
    % tag blobs with dates:
    timetag = cellstr(datestr(blob_time,'mmmdd'));
    %text(SSTa, AxisRatio, timetag,'fontsize',11);
    hb = colorbar;
    set(get(hb,'xlabel'),'string','background SST (dgr)');
    xlabel('SST anomaly (K)');
    ylabel('Major to Minor Axis Ratio');
    
end
set(gca,'fontsize',14);
grid on
title({'compare blob characters on days','with and without SWC'});

% It seems that the days without SWC tends to have smaller blobs with
% weaker SST anoamly. 

% group the RHB samples by day as well (with and without significant
% coherence)

%% 3. make distribution plots
% this one is what I need.
figure(3); clf;
nval=4;
spatialres=2;   %km
for j = 1:2
    flag = flagNames{j};
    
    %stats_grp = get_blob_characters(RHBblobs.(flag), spatialres);
    %varlist2 = fieldnames(stats_grp);
    stats_grp.Strength = [RHBblobs.(flag).max_SSTa];
    tmp_stats =  [RHBblobs.(flag).stats_selected];
    stats_grp.Shape = unique([tmp_stats.MajAxisLen])./[tmp_stats.MinAxisLen];
    stats_grp.Size = [tmp_stats.EqvDiam]*spatialres;                                   % what is the units?
    stats_grp.backgroundSST = [RHBblobs.(flag).ave_SSTbg]-273.15;
    %SSTbg = [RHBblobs.(flag).ave_SSTbg]-273.15;
    %blob_time = [RHBblobs.(flag).time];
    
    for iv = 1:nval
        VN = varlist{iv};
        val = stats_grp.(VN);
        
        ip = iv+nval*(j-1);
        subplot(2,nval,ip);
        histogram(val, edges.(VN));
        hold on
        yrange = get(gca,'ylim');
        plot(median(val)*ones(1,2),yrange,'color', RHBcolor,'linewidth',1.2);
        if iv ==2
            set(gca,'XTick',[0:20:200]);
        end
        xlabel(xlabels{iv});
        ylabel('count');
        title(VN);
    end
    
end
xc_savefig(gcf,'Figs','Compare_blob_key_characters_for_days_with_SWC_and_withoutSWC_addSSTbg.jpg',[0 0 12 8]);






