% explore:
% % try: finding 2D histogram of daily cloudiness and u_vec * SST grad (change of SST experienced by an air parcel).
% Attempt1 : u_vec: ERA5 daily averages
%            SST gradients: L4 daily averaged.

clear all; clc;

addpath('/Users/xchen/Documents/GitHub/SAM_LES/matlab');
% read the following:
dataroot = pwd;     % make sure to execute the script from the right path.

dataFolder.GOES = 'gpblend'; dataFolder.ERA5 = 'era5_data';
ERA5_dataFN = 'era5_10mWindVectors_DJF2019-2022.nc';
L4_dataFN = 'gpblend_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N.mat';

load([dataroot filesep dataFolder.GOES filesep L4_dataFN]);
L4data = data;  clear data

absFN = [dataroot filesep dataFolder.ERA5 filesep ERA5_dataFN];
ERA5data=read_netCDF_into_matlab_structure(absFN);
ERA5data.wdir = atan2(ERA5data.vwnd, ERA5data.uwnd);


% compute SST gradient from L4 data: 
cutoff_scale_km = 600;
for it = 1:length(L4data.time)
    SST0 = L4data.analysed_sst(:,:,it)';
    lon = double(L4data.lon(:,it));
    lat = double(L4data.lat(:,it));
    
    midlat = (lat(1:end-1) + lat(2:end))/2;
    
    xvec = [0; cumsum(diff(lon)*111E3.*cosd(midlat))];   % km
    yvec = [0; cumsum(diff(lat)*111E3)];                 % km
    
    meanSST_LS(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',false);    % the size is different here.
    
    SST_anom_spatial = SST0 - meanSST_LS;
%     SSTgrad(it) = spatial_grad(SST_anom_spatial, xvec, yvec);  %units: K/km
%     SSTLaplacian(:,:, it) = spatial_laplacian(SST_anom_spatial, xvec, yvec);
    
    SSTgrad(it) = spatial_grad(SST0, xvec, yvec);  %units: K/km
    SSTLaplacian(:,:, it) = spatial_laplacian(SST0, xvec, yvec);
    
end

[LON_L4, LAT_L4] = meshgrid(L4data.lon(:,1), L4data.lat(:,1));
LON = LON_L4; LAT = LAT_L4;
% select subregion:
subregion = [-59, -48; 8 18];
lonmask = L4data.lon(:,1)>=subregion(1,1) & L4data.lon(:,1)<=subregion(1,2);
latmask = L4data.lat(:,1)>=subregion(2,1) & L4data.lat(:,1)<=subregion(2,2);

LON_sub = LON_L4(latmask, lonmask); 
LAT_sub = LAT_L4(latmask, lonmask); 


cmap = getPyPlot_cMap('RdBu_r',64);
% visualize gradient magnitude:
figure(1)
for it = 1:length(L4data.time)
    
    %pcolor(lon, lat, SSTLaplacian(:,:,it)); shading flat;
    pcolor(lon, lat, SSTgrad(it).mag); shading flat;
    
    colorbar
    colormap(parula);
    hold on
    quiver(LON(1:4:end, 1:4:end),LAT(1:4:end, 1:4:end), ...
        SSTgrad(it).xcomp(1:4:end, 1:4:end), SSTgrad(it).ycomp(1:4:end, 1:4:end), 1, 'k');
    hold off
    
    title(datestr(L4data.time_num(it)));
    pause(0.4);
    
end

% compute uvec * SST grad 
% deal with resolution:
[ERA5data.LON, ERA5data.LAT] = meshgrid(ERA5data.lon, ERA5data.lat);
dx = diff(LON,1, 2);
dx = [zeros(size(dx,1),1)  dx];

dy = diff(LAT,1,1);
dy = [zeros(1,size(dy,2)); dy];
XX = 111*cosd(LAT).*cumsum(dx,2);     % km
YY = 111.*cumsum(dy,1);                        % km
for it = 1:length(ERA5data.time)
    ERA5_uc(:,:,it) = interp2(ERA5data.LON-360, ERA5data.LAT, ERA5data.uwnd(:,:,it)', LON, LAT);
    ERA5_vc(:,:,it) = interp2(ERA5data.LON-360, ERA5data.LAT, ERA5data.vwnd(:,:,it)', LON, LAT);
    ERA5_wdiv(:,:,it) = divergence(XX, YY, ERA5_uc(:,:,it), ERA5_vc(:,:,it));
end


for it = 1:length(L4data.time_num)
    tid_era5 = find(ERA5data.time_num==floor(L4data.time_num(it)));
    air_Ttend(:,:,it) = (ERA5_uc(:,:,tid_era5).*SSTgrad(it).xcomp + ERA5_vc(:,:,tid_era5) .* SSTgrad(it).ycomp)*3.6;  % units: K/h; (1m/s = 3.6km/h)
    winddiv(:,:,it) = ERA5_wdiv(:,:,tid_era5);
end

%% daily cloudiness:
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
    
    SST_all = permute(GOES_ATOMIC.sea_surface_temperature, [2,1,3]);
    time_all = GOES_ATOMIC.time_num;
    
    % check cloud mask and throw away maps that contains all cloudy pixels.
    %% QC data:
    % through away maps that is 100% garbage;
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
    
end

% cldfreq is on a finer 2km grid: coarsen it up:
for it = 1:length(L4data.time)
    cldfreq_coarse(:,:,it) = interp2(GOES_LON, GOES_LAT, cldfreq(:,:,it), LON, LAT,'natural');
end


% how to plot x as percentiles?
% 1. create edges based on percentiles:
mean_air_Ttend = mean(air_Ttend(latmask, lonmask,:), 3, 'omitnan');
mean_cldfreq = mean(cldfreq_coarse(latmask, lonmask,:),3,'omitnan');

prcs = [5:5:95];
aTedges = prctile(mean_air_Ttend(:), prcs);
aTedges_wider = prctile(mean_air_Ttend(:), [10:10:90]);
CFedges = prctile(mean_cldfreq(:), [10:10:90]);

% now: plot joint distribution;
% aTedges = [-8:0.1:8].*10^-4;
% aTedges_wider = [-3:1:3].*10^(-4);
% CFedges = [0:0.1:1];

figure(3); clf;
hall = histogram2(mean_air_Ttend, mean_cldfreq, ...
    aTedges, CFedges, 'DisplayStyle','tile','ShowEmptyBins','on');
colorbar;
occur_counts_all = hall.Values;
hold on
plot(mean(mean_air_Ttend(:)), mean(mean_cldfreq(:)),'+r','markersize',12,'linewidth',2);
hold on
plot(median(mean_air_Ttend(:)), median(mean_cldfreq(:)), 'dr','markersize',12,'linewidth',2);


% now look at the distribution according to different large scale
% conditions:
dayType = get_dayType_by_U10_and_LTS;
cloudtypes = fieldnames(dayType);
for tt = 1:4
    CN = cloudtypes{tt};
    tids = ismember(floor(L4data.time_num), dayType.(CN));
    aTtend.(CN) = air_Ttend(latmask, lonmask, tids);
    CFc.(CN) = cldfreq_coarse(latmask, lonmask, tids);
    SSTlap.(CN) = SSTLaplacian(latmask, lonmask, tids);
    wdiv.(CN) = winddiv(latmask, lonmask, tids);
    SST_LS.(CN) = meanSST_LS(latmask,lonmask, tids);
    
end


% plot effective SST gradient into a movie:
hfig = figure(1); clf;
for tt = 1:4
    CN = cloudtypes{tt};
    clear images
    for it = 1:size(aTtend.(CN),3)
        
        pcolor(LON_sub, LAT_sub, aTtend.(CN)(:,:,it)); shading flat;
        
        colorbar
        colormap(redblue);
        hold on
        contour(LON_sub, LAT_sub, aTtend.(CN)(:,:,it),'k');
        %         quiver(LON(1:4:end, 1:4:end),LAT(1:4:end, 1:4:end), ...
        %             SSTgrad(it).xcomp(1:4:end, 1:4:end), SSTgrad(it).ycomp(1:4:end, 1:4:end), 1, 'k');
        %         hold off
        
        title(datestr(L4data.time_num(it)));
        caxis([-8*1e-4, 8*1e-4]);
        pause(0.1);
        frame = getframe(hfig);
        images{it} = frame2im(frame);
        
    end
    gifname = [CN '_effective_SSTgrad.gif'];
    gif_abspath = ['Figs' filesep gifname];
    write_frames_into_video(images, gif_abspath, 'gif');
    
end


figure(4); clf;
for tt = 1:4
    CN = cloudtypes{tt};
    subplot(2,2,tt)
    
    % find the percentile first:
    aTedges = prctile(mean(aTtend.(CN),3, 'omitnan'), [prcs], 'all');
    aTedges_wider = prctile(mean(aTtend.(CN),3, 'omitnan'), [0:10:100],'all');
    CFedges = prctile(mean(CFc.(CN),3,'omitnan'), [10:10:90],'all');

    
    hh(tt) = histogram2(mean(aTtend.(CN),3, 'omitnan'), mean(CFc.(CN),3,'omitnan'), ...
        aTedges_wider, CFedges, 'DisplayStyle','tile','ShowEmptyBins','off','normalization','count');
%     
%     hh = histogram2(aTtend.(CN), CFc.(CN), ...
%         aTedges, CFedges, 'DisplayStyle','tile','ShowEmptyBins','off','normalization','count');

    colorbar;
    %caxis([0 0.05])
    hold on
    plot(mean(mean(aTtend.(CN),3, 'omitnan'), 'all'), mean(mean(CFc.(CN),3,'omitnan'),'all'),'+r','markersize',12,'linewidth',2);

    plot(median(mean(aTtend.(CN),3, 'omitnan'),'all'), median(mean(CFc.(CN),3,'omitnan'),'all'), 'dr','markersize',12,'linewidth',2);
    %caxis([0 350]);
    
    %occur_counts{tt}=hh.Values;
    xlabel('\bf{u} \cdot \nablaSST (K/hr)');
    ylabel('cloudiness');
    set(gca,'fontsize',13);
end

SSTLpEdges=[-0.2:0.02:0.2]*10^(-8);
SSTLpEdges_wider =[-1.5:0.5:1.5].*10^(-9);

figure(5)
for tt = 1:4
    CN = cloudtypes{tt};
    
    SSTLpEdges = prctile(mean(SSTlap.(CN),3,'omitnan'), prcs, 'all');
    SSTLpEdges_wider = prctile(mean(SSTlap.(CN),3,'omitnan'), [10:10:90],'all');
    
    subplot(2,2,tt)
    hh_SSTlap(tt) = histogram2(mean(SSTlap.(CN),3,'omitnan'), mean(CFc.(CN),3,'omitnan'), ...
        SSTLpEdges, CFedges, 'DisplayStyle','tile','ShowEmptyBins','on','normalization','count');
    colorbar;
   % caxis([0 0.05])
    
    %occur_counts_laplacian{tt}=hh.Values;
    
end


% check if the ensemble averaged aTend and cldfreq changes:
pids = [2 3 5 6];
figure(6)
subplot(2,3,[1 4]);
pcolor(LON, LAT, mean(air_Ttend, 3,'omitnan')); shading flat;
colorbar;
colormap(cmap)
axis('square');
xlim([-59 -48]); ylim([8 18]);
caxis([-4 4]*10^-4)

for ip = 1:4
    subplot(2,3,pids(ip))
    CN = cloudtypes{ip};
    pcolor(LON(latmask, lonmask), LAT(latmask, lonmask), mean(aTtend.(CN), 3,'omitnan')); shading flat;
    xlim([-59 -48]); ylim([8 18]);
    colorbar;
    axis('square');
    colormap(cmap)
    caxis([-4 4]*10^-4)

end

%  find out where we see increased percentage of negative SST laplacian for
%  sugar condition:
CN='sugar';
selmask = mean(SSTlap.(CN),3,'omitnan')<=-0.5*10^(-9) & mean(CFc.(CN),3,'omitnan')>=0.4;  % base on the normalized occurrence. 
LONs = repmat(LON(latmask, lonmask),1,1);
LATs = repmat(LAT(latmask, lonmask),1,1);

dataIn = SSTlap.(CN);
levs = [-2.5:0.5:2.5].*10^-9;
show_selected_locations(LONs, LATs, dataIn, selmask,levs);

tail_SSTlap_locs.lon = LONs(selmask);
tail_SSTlap_locs.lat = LATs(selmask);

figure
pcolor(LONs, LATs,mean(SSTlap.(CN),3,'omitnan')*10^2); shading flat;
colorbar; caxis([-1.6, 1.6]*10^(-7))
hold on
plot(tail_SSTlap_locs.lon, tail_SSTlap_locs.lat,'.m','linewidth',0.2,'markersize',0.1);
title([CN '-favored day averaged'])

%histogram2(tail_SSTlap_locs.lon, tail_SSTlap_locs.lat, [-59:0.25:-48],[8:0.25:18], 'DisplayStyle','tile','ShowEmptyBins','off','normalization','count');
%colorbar;

%% air rate of chagne of SST
%aTedges_wider = [-3:1:3]*10^(-4);
%CFedges_narrower = [0:0.05:1];
% Xedges = aTedges;
% Yedges = CFedges; 

Xedges = [0:10:100];
Yedges = [0:10:100];

lncols_pr = getPyPlot_cMap('jet',length(Yedges)-1);

pos = customize_subplot_size(2,3,0.1,0.1);
figure
hsub0 = subplot(2,3,[1,4]);
hax1 = gca;
Y.val = mean(cldfreq_coarse(latmask,lonmask,:), 3,'omitnan');
X.val = mean(air_Ttend(latmask,lonmask,:),3,'omitnan') ;
Y.labelstr = 'cloudiness';
X.labelstr = '\bf{u} \cdot \nablaSST (K/hr)';

yscale = 'log';

hlgd = plot_probability_results(hax1, X, Y, Xedges, Yedges, 'percentile', lncols_pr, yscale);
%xlim([-4,4]*10^-4)
axis('square')
set(hlgd,'loc','northwest');
%ylim(10.^[-5,0]);



for tt =1:4
    CN = cloudtypes{tt};
    
    hsub(tt) =subplot(2,3,pids(tt));
    hax = gca;
    Y.val = mean(CFc.(CN), 3,'omitnan');
    X.val = mean(aTtend.(CN),3,'omitnan') ;
    
    hlgd = plot_probability_results(hax, X, Y, Xedges, Yedges, 'percentile', lncols_pr, yscale);
    
    set(hsub(tt), 'position',pos{pids(tt)});
    set(hlgd,'loc','northwest');
    %ylim([0 0.45]);
    %xlim([-4,4]*10^-4)
    %ylim(10.^[-5,0]);
    title([CN '-favored day averaged'])
end


%% SST laplacian:
Xedges = [0:10:100]; %SSTLpEdges*10^2;
Yedges = [0:10:100]; %CFedges; 

lncols_pr = getPyPlot_cMap('jet',length(Yedges)-1);

pos = customize_subplot_size(2,3,0.1,0.1);
figure
hsub0 = subplot(2,3,[1,4]);
hax1 = gca;
Y.val = mean(cldfreq_coarse(latmask,lonmask,:), 3,'omitnan');
X.val = mean(SSTLaplacian(latmask,lonmask,:),3,'omitnan')*10^2 ;
Y.labelstr = 'cloudiness';
X.labelstr =  '{\nabla}^2SST (K/(10km)^2)';


yscale = 'log';

hlgd = plot_probability_results(hax1, X, Y, Xedges, Yedges, 'percentile', lncols_pr, yscale);
%xlim([-1.6 1.6]*10^(-9+2));
axis('square')
set(hlgd,'loc','northwest');
%ylim(10.^[-5,0]);



for tt =1:4
    CN = cloudtypes{tt};
    
    hsub(tt) =subplot(2,3,pids(tt));
    hax = gca;
    Y.val = mean(CFc.(CN), 3,'omitnan');
    X.val = mean(SSTlap.(CN),3,'omitnan')*10^2;
    
    hlgd = plot_probability_results(hax, X, Y, Xedges, Yedges, 'percentile', lncols_pr, yscale);
    
    set(hsub(tt), 'position',pos{pids(tt)});
    set(hlgd,'loc','northwest');
    %ylim([0 0.45]);
    %xlim([-4,4]*10^-4)
   % ylim(10.^[-5,0]);
    %xlim([-1.6 1.6]*10^(-9+2));

    title([CN '-favored day averaged'])
end





%% find out where swe see the tail values of rate of change in SST for sugar and gravel condition.
% compute standard deviation:
CN='gravel';
selmask = (mean(aTtend.(CN),3, 'omitnan'))>=1.5*10^(-4) & mean(CFc.(CN),3,'omitnan')>=0.4;  % base on the normalized occurrence. 
LONs = repmat(LON(latmask, lonmask),1,1);
LATs = repmat(LAT(latmask, lonmask),1,1);

tail_aTend_locs.lon = LONs(selmask);
tail_aTend_locs.lat = LATs(selmask);

figure;
subplot(1,2,1);
pcolor(LONs, LATs,mean(aTtend.(CN),3, 'omitnan')); shading flat;
hold on
contour(LONs, LATs, std(aTtend.gravel,1,3),[1:0.5:2.5].*10^(-4),'k');
hb =colorbar; caxis([-2.5, 2.5]*10^(-4))
set(get(hb,'xlabel'), 'string', '\bf{u} \cdot \nablaSST (K/hr)');
plot(tail_aTend_locs.lon, tail_aTend_locs.lat,'.m', 'linewidth',0.2,'markersize',0.1);
title([CN '-favored days averaged'])
xlabel('longitude'); ylabel('latitude');
set(gca,'fontsize',14);
axis('equal')

subplot(1,2,2);
pcolor(LONs, LATs,std(aTtend.gravel,1,3)); shading flat;
hold on
contour(LONs, LATs, std(aTtend.gravel,1,3),[1:0.5:2.5].*10^(-4),'k');
hb =colorbar; caxis([0, 2.5]*10^(-4));
set(get(hb,'xlabel'), 'string', '\bf{u} \cdot \nablaSST (K/hr)');

plot(tail_aTend_locs.lon, tail_aTend_locs.lat,'.m', 'linewidth',0.2,'markersize',0.1);
title(['standard deviation'])
xlabel('longitude'); ylabel('latitude');
set(gca,'fontsize',14);
axis('equal')


%histogram2(tail_aTend_locs.lon, tail_aTend_locs.lat, [-59:0.25:-48],[8:0.25:18], 'DisplayStyle','tile','ShowEmptyBins','off','normalization','count');
%colorbar;


CN='sugar';
selmask = (mean(aTtend.(CN),3, 'omitnan'))>1.5*10^(-4) & mean(CFc.(CN),3,'omitnan')>=0.4;  % base on the normalized occurrence. 
selmask2 = mean(aTtend.(CN),3, 'omitnan')<-0.5*10^(-4) & mean(CFc.(CN),3,'omitnan')>=0.7;
selmask3 = mean(aTtend.(CN),3, 'omitnan')>0.5*10^(-4) & mean(CFc.(CN),3,'omitnan')>=0.7;

tail_aTend_locs.lon = LONs(selmask);
tail_aTend_locs.lat = LATs(selmask);


figure
subplot(1,2,1);
pcolor(LONs, LATs,mean(aTtend.(CN),3, 'omitnan')); shading flat;
colorbar; caxis([-4, 4]*10^(-4))
hold on
hb =colorbar; caxis([-2.5, 2.5]*10^(-4))
set(get(hb,'xlabel'), 'string', '\bf{u} \cdot \nablaSST (K/hr)');

plot(tail_aTend_locs.lon, tail_aTend_locs.lat,'.m', 'linewidth',0.2);
contour(LONs, LATs, std(aTtend.(CN),1,3),[1:0.5:2.5].*10^(-4),'k');

%plot(LONs(selmask2), LAT(selmask2),'.c');
%plot(LONs(selmask3), LAT(selmask3),'.y');
title([CN '-favored day averaged'])

xlabel('longitude'); ylabel('latitude');
set(gca,'fontsize',14);
axis('equal')

subplot(1,2,2);
pcolor(LONs, LATs,std(aTtend.(CN),1,3)); shading flat;
hold on
plot(tail_aTend_locs.lon, tail_aTend_locs.lat,'.m', 'linewidth',0.2,'markersize',0.01);

contour(LONs, LATs, std(aTtend.(CN),1,3),[1:0.5:2.5].*10^(-4),'k');
hb =colorbar; caxis([0, 2.5]*10^(-4));
set(get(hb,'xlabel'), 'string', '\bf{u} \cdot \nablaSST (K/hr)');

title(['standard deviation'])
xlabel('longitude'); ylabel('latitude');
set(gca,'fontsize',14);
axis('equal')



histogram2(tail_aTend_locs.lon, tail_aTend_locs.lat, [-59:0.25:-48],[8:0.25:18], 'DisplayStyle','tile','ShowEmptyBins','off','normalization','count');
colorbar;


%% figure: plot averaged u*SSTgrad for the gravel case, and highlight location where there are extreme values:



%%
% find probability for cloud bins conditioned on different aTtend or SST Laplacian:
aTedges_wider =  [0:25:100]; %[-3:1:3]*10^(-4);
CFedges_narrower = [0:10:100]; %[0:0.05:1];
lncols_rb = getPyPlot_cMap('RdBu_r',length(aTedges_wider)-1);

pos = customize_subplot_size(2,3,0.1,0.1);
figure(22)
hsub0 = subplot(2,3,[1,4]);
hax1 = gca;
X.val = mean(cldfreq_coarse(latmask,lonmask,:), 3,'omitnan');
Y.val = mean(air_Ttend(latmask,lonmask,:),3,'omitnan') ;
X.labelstr = 'cloudiness';
Y.labelstr = '\bf{u} \cdot \nablaSST (K/hr)';

yscale = 'linear';

Xedges = CFedges_narrower; 
Yedges = aTedges_wider;
plot_probability_results(hax1, X, Y, Xedges, Yedges, 'percentile', lncols_rb, yscale);


for tt =1:4
    CN = cloudtypes{tt};
    
    hsub(tt) =subplot(2,3,pids(tt));
    hax = gca;
    X.val = mean(CFc.(CN), 3,'omitnan');
    Y.val = mean(aTtend.(CN),3,'omitnan') ;
    
    plot_probability_results(hax, X, Y, Xedges, Yedges, 'percentile',lncols_rb, yscale);
    
    %set(hsub(tt), 'position',pos{pids(tt)});
    %ylim([0 0.45]);
    title([CN '-favored day averaged'])
end



% use the new code to find probability for SST laplacian:
CFedges_narrower =  [0:10:100]; %[0:0.1:1];
SSTLpEdges_wider =  [0:20:100]; %[-1.5:0.25:1.5].*10^(-9);

Xedges = CFedges_narrower; 
Yedges = SSTLpEdges_wider;

lncols_rb = getPyPlot_cMap('coolwarm',(length(Yedges)-1)*1);
mid = round(length(lncols_rb)/2);
lncols_rb(mid,:) = [0.45, 0.45, 0.45];


pos = customize_subplot_size(2,3,0.1,0.1);
figure
hsub0 = subplot(2,3,[1,4]);
hax1 = gca;
X.val = mean(cldfreq_coarse(latmask,lonmask,:), 3,'omitnan');
Y.val = mean(SSTLaplacian(latmask,lonmask,:),3,'omitnan') ;
X.labelstr = 'cloudiness';
Y.labelstr = '\nabla^2SST (K/(10km)^2)';

yscale = 'linear';

plot_probability_results(hax1, X, Y, Xedges, Yedges, 'percentile',lncols_rb, yscale);

for tt =1:4
    CN = cloudtypes{tt};
    
    hsub(tt) =subplot(2,3,pids(tt));
    hax = gca;
    X.val = mean(CFc.(CN), 3,'omitnan');
    Y.val = mean(SSTlap.(CN),3,'omitnan') ;
    
    plot_probability_results(hax, X, Y, Xedges, Yedges,'percentile', lncols_rb, yscale);
    
    set(hsub(tt), 'position',pos{pids(tt)});
    %ylim([0 0.45]);
    title([CN '-favored day averaged'])
end


%% test do no averaging: a lot harder to see information (even though the magnitude of xvar will be larger.
Xedges = aTedges;
Yedges = CFedges(1:2:end); 

lncols_pr = getPyPlot_cMap('Paired',length(Yedges)-1);

pos = customize_subplot_size(2,3,0.1,0.1);
figure
hsub0 = subplot(2,3,[1,4]);
hax1 = gca;
Y.val = cldfreq_coarse(latmask,lonmask,:);
X.val = air_Ttend(latmask,lonmask,:) ;
Y.labelstr = 'cloudiness';
X.labelstr = '\bf{u} \cdot \nablaSST (K/hr)';

yscale = 'log';

hlgd = plot_probability_results(hax1, X, Y, Xedges, Yedges, lncols_pr, yscale);
%xlim([-4,4]*10^-4)
axis('square')
set(hlgd,'loc','northwest');
%ylim(10.^[-5,0]);



for tt =1:4
    CN = cloudtypes{tt};
    
    hsub(tt) =subplot(2,3,pids(tt));
    hax = gca;
    Y.val = CFc.(CN);
    X.val = aTtend.(CN) ;
    
    hlgd = plot_probability_results(hax, X, Y, Xedges, Yedges, lncols_pr, yscale);
    
    set(hsub(tt), 'position',pos{pids(tt)});
    set(hlgd,'loc','northwest');
    %ylim([0 0.45]);
    %xlim([-4,4]*10^-4)
   % ylim(10.^[-5,0]);
    title([CN '-favored day averaged'])
end



%% test function and the relative change in a similar way as in Desbin et al. (2021)
figure(10); clf
for tt = 1:4
    CN = cloudtypes{tt};
    [RCC_l, RCC_dm] = compute_relative_change_of_cloudiness_in_xbins(aTtend.(CN), [0:5:100], CFc.(CN), mean_cldfreq, LON_sub, LAT_sub);
    null(tt) = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(aTtend.(CN), 5, CFc.(CN), mean_cldfreq, 5000);

    xbin = 2.5:5:100;
     %[RCC_l, RCC_dm] = compute_relative_change_of_cloudiness_in_xbins(mean(aTtend.(CN),3), [0:5:100], mean(CFc.(CN),3,'omitnan'), mean_cldfreq, LON_sub, LAT_sub);

    %null(tt) = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(mean(aTtend.(CN),3), 5, mean(CFc.(CN),3,'omitnan'), mean_cldfreq, 1000);

    
   % prcs_1 = prctile(mean(aTtend.(CN),3), [0:5:100],'all');
    prcs_1 = prctile(aTtend.(CN), [0:5:100],'all');
    prc_rank = interp1(prcs_1, [0:5:100], 0);
    subplot(2,2,tt);
    hold on
    %bar(xbin, RCC_dm, 'FaceAlpha',1);
    bar(xbin, RCC_l, 1.0,'FaceAlpha', 0.5);
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
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_cloudiness_vs_effective_SSTgrad_percentiles_4types_v2.jpg', [0 0 10 8]); 

% test: 
selcrit = aTtend.(CN)<= prctile(aTtend.(CN), 10);
sel_wdiv = wdiv.(CN)(selcrit);
wconv_frac = numel(find(sel_wdiv<0))/numel(sel_wdiv);

%% effective SST gradient:
[RCC_l, RCC_dm, xedge] = compute_relative_change_of_cloudiness_in_xbins(air_Ttend(latmask, lonmask,:), [0:5:100], cldfreq_coarse(latmask, lonmask,:), mean_cldfreq, LON_sub, LAT_sub);
xbin = 2.5:5:100;

% do an estimation of 5% null hypothesis level:
% randomly took 5% of data 1000 times, and find the standard deviation of
% average cloudiness; use 1.96*stdv as the range.
null_results = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(air_Ttend(latmask, lonmask,:), 5, cldfreq_coarse(latmask, lonmask,:), mean_cldfreq);

    
prcs_1 = prctile(air_Ttend(latmask, lonmask,:), [0:5:100],'all');
prc_rank = interp1(prcs_1, [0:5:100], 0);
figure(11); clf;
hold on
%bar(xbin, RCC_dm, 'FaceAlpha',1);
bar(xbin, RCC_l,1.0, 'FaceAlpha', 0.5);
plot([prc_rank, prc_rank],[-15, 30], '--b', 'linewidth', 1.2);
hold on;
% plot 95% confidence level of the mean change of CF in 5% of data; 
vlocs = [0, null_results.stdv*1.96 + null_results.mean; 
         100, null_results.stdv*1.96 + null_results.mean; 
         100, -null_results.stdv*1.96 + null_results.mean; 
         0, -null_results.stdv*1.96 + null_results.mean; ];
         
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
    
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_cloudiness_vs_effective_SSTgrad_percentiles.jpg', [0 0 10 8]); 


%% wind div:
[RCC_l, RCC_dm, xedge] = compute_relative_change_of_cloudiness_in_xbins(winddiv(latmask, lonmask,:), [0:5:100], cldfreq_coarse(latmask, lonmask,:), mean_cldfreq, LON_sub, LAT_sub);
xbin = 2.5:5:100;

% do an estimation of 5% null hypothesis level:
% randomly took 5% of data 1000 times, and find the standard deviation of
% average cloudiness; use 1.96*stdv as the range.
null_results = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(winddiv(latmask, lonmask,:), 5, cldfreq_coarse(latmask, lonmask,:), mean_cldfreq);

    
prcs_1 = prctile(winddiv(latmask, lonmask,:), [0:5:100],'all');
prc_rank = interp1(prcs_1, [0:5:100], 0);
figure(1); clf;
hold on
%bar(xbin, RCC_dm, 'FaceAlpha',1);
bar(xbin, RCC_l,1.0, 'FaceAlpha', 0.5);
hold on;
% plot 95% confidence level of the mean change of CF in 5% of data; 
vlocs = [0, null_results.stdv*1.96 + null_results.mean; 
         100, null_results.stdv*1.96 + null_results.mean; 
         100, -null_results.stdv*1.96 + null_results.mean; 
         0, -null_results.stdv*1.96 + null_results.mean; ];
         
patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');    
yrange = get(gca, 'ylim');
plot([prc_rank, prc_rank],yrange, '--b', 'linewidth', 1.2);

xlim([0 100]);    
xlabel('\nabla \cdot \bf{u} percentiles');
ylabel('average change of cloudiness in each percentiles (%)');
set(gca,'fontsize',14);
title('Jan and Feb 2020 (all)');
grid on
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_cloudiness_vs_winddiv_percentiles.jpg', [0 0 10 8]); 


figure(2); clf
for tt = 1:4
    CN = cloudtypes{tt};
    [RCC_l, RCC_dm] = compute_relative_change_of_cloudiness_in_xbins(wdiv.(CN), [0:5:100], CFc.(CN), mean_cldfreq, LON_sub, LAT_sub);
    xbin = 2.5:5:100;
    
    null(tt) = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(wdiv.(CN), 5, CFc.(CN), mean_cldfreq, 5000);
    
    
    prcs_1 = prctile(wdiv.(CN)(:), [0:5:100]);
    prc_rank = interp1(prcs_1, [0:5:100], 0);
    subplot(2,2,tt);
    hold on
    %bar(xbin, RCC_dm, 'FaceAlpha',1);
    bar(xbin, RCC_l, 'FaceAlpha', 0.5);
    plot([prc_rank, prc_rank],[-10, 15], '--k');
    hold on;
    % plot 95% confidence level of the mean change of CF in 5% of data;
    vlocs = [0, null(tt).stdv*1.96 + null(tt).mean;
        100, null(tt).stdv*1.96 + null(tt).mean;
        100, -null(tt).stdv*1.96 + null(tt).mean;
        0, -null(tt).stdv*1.96 + null(tt).mean; ];
    
    patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');
    xlim([0 100]);
    xlabel('\nabla \cdot \bf{u} percentiles');
    ylabel('average change of cloudiness in each percentiles (%)');
    set(gca,'fontsize',14);
    title(CN);
    grid on
end
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_cloudiness_vs_effective_SSTgrad_percentiles_4types.jpg', [0 0 10 8]); 


%% SST laplacian:
[RCC_l, RCC_dm, xedge2] = compute_relative_change_of_cloudiness_in_xbins(SSTLaplacian(latmask, lonmask,:), [0:5:100], cldfreq_coarse(latmask, lonmask,:), mean_cldfreq, LON_sub, LAT_sub);
xbin = 2.5:5:100;

% do an estimation of 5% null hypothesis level:
% randomly took 5% of data 1000 times, and find the standard deviation of
% average cloudiness; use 1.96*stdv as the range.
null_lap = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(SSTLaplacian(latmask, lonmask,:), 5, cldfreq_coarse(latmask, lonmask,:), mean_cldfreq);

    
prcs_1 = prctile(SSTLaplacian(latmask, lonmask,:), [0:5:100],'all');
prc_rank = interp1(prcs_1, [0:5:100], 0);
figure(12); clf;
hold on
%bar(xbin, RCC_dm, 'FaceAlpha',1);
bar(xbin, RCC_l,1.0, 'FaceAlpha', 0.5);
plot([prc_rank, prc_rank],[-15, 30], '--b', 'linewidth', 1.2);
hold on;
% plot 95% confidence level of the mean change of CF in 5% of data; 
vlocs = [0, null_lap.stdv*1.96 + null_lap.mean; 
         100, null_lap.stdv*1.96 + null_lap.mean; 
         100, -null_lap.stdv*1.96 + null_lap.mean; 
         0, -null_lap.stdv*1.96 + null_results.mean; ];
         
patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');    
xlim([0 100]);    
xlabel('{\nabla}^2SST percentiles');
ylabel('average change of cloudiness in each percentiles (%)');
set(gca,'fontsize',14);
title('Jan and Feb 2020 (all)');
grid on
xc_savefig(gcf, './Figs/probability_exploration', 'change_of_cloudiness_vs_SSTLaplacian_percentiles.jpg', [0 0 10 8]); 


figure(13); clf
for tt = 1:4
    CN = cloudtypes{tt};
    [RCC_l, RCC_dm] = compute_relative_change_of_cloudiness_in_xbins(SSTlap.(CN), [0:5:100], CFc.(CN), mean_cldfreq, LON_sub, LAT_sub);
    xbin = 2.5:5:100;
    
    null(tt) = bootstrap_nullhypothesis_level_for_relative_cloudiness_change(SSTlap.(CN), 5, CFc.(CN), mean_cldfreq, 1000);

    
    prcs_1 = prctile(SSTlap.(CN)(:), [0:5:100]);
    prc_rank = interp1(prcs_1, [0:5:100], 0);
    subplot(2,2,tt);
    hold on
    %bar(xbin, RCC_dm, 'FaceAlpha',1);
    bar(xbin, RCC_l, 'FaceAlpha', 0.5);
    plot([prc_rank, prc_rank],[-10, 15], '--k');
    hold on;
    % plot 95% confidence level of the mean change of CF in 5% of data;
    vlocs = [0, null(tt).stdv*1.96 + null(tt).mean;
        100, null(tt).stdv*1.96 + null(tt).mean;
        100, -null(tt).stdv*1.96 + null(tt).mean;
        0, -null(tt).stdv*1.96 + null(tt).mean; ];
    
    patch(vlocs(:,1), vlocs(:,2), [0.55, 0.55, 0.55],'FaceAlpha', 0.3,'EdgeColor','r');
    xlim([0 100]);
    xlabel('{\nabla}^2SST percentiles');
    ylabel('average change of cloudiness (%)');
    set(gca,'fontsize',14);
    title(CN);
    grid on
end

%% save the locations
% find locations of SST grad that is at the end of the distribution:

% I guess I can show this.
CN='fish';
clear wanom_pdf canom_pdf pdf_levs
for tt = 1:4
    CN = cloudtypes{tt};
    CFRC.(CN) = (CFc.(CN) - mean_cldfreq)./mean_cldfreq;
    crit1 = aTtend.(CN)<= prctile(aTtend.(CN)(:),5);
    CFRC_thres1 = prctile(CFRC.(CN)(crit1),75,'all');
    
    crit2 = aTtend.(CN)> prctile(aTtend.(CN)(:),95);
    CFRC_thres2 = prctile(CFRC.(CN)(crit2),75,'all');
    
    % i want to select the extreme gradient that is less variable in time.
    varcrit3 = abs(aTtend.(CN)-mean(aTtend.(CN),3,'omitnan') )<= 1.0*std(aTtend.(CN), 1, 3,'omitnan');
    %varcrit4 = abs(aTtend.(CN)-mean(aTtend.(CN),3,'omitnan') )> 3*std(aTtend.(CN), 1, 3,'omitnan');
    
    % CFRC: in the upper quartiles:
    selmask = aTtend.(CN)<= prctile(aTtend.(CN)(:),5) & varcrit3;   %& CFRC.(CN)>=CFRC_thres1; % & wdiv.(CN)<=prctile(wdiv.(CN)(:),5);  % base on the normalized occurrence.
    selmask2 = aTtend.(CN)> prctile(aTtend.(CN)(:),95) & varcrit3; %& CFRC.(CN)>=CFRC_thres2; % & wdiv.(CN)<=prctile(wdiv.(CN)(:),5);   % base on the normalized occurrence.
    
    LONs = repmat(LON_sub, 1,1, size(aTtend.(CN),3));
    LATs = repmat(LAT_sub, 1,1, size(LONs,3));
    
    tail_aTend_locs.lon = LONs(selmask);
    tail_aTend_locs.lat = LATs(selmask);
    
    lon_edges = [-59:0.25:-48]; lat_edges=[8:0.25:18];
    lon_bin = 0.5*(lon_edges(1:end-1)+lon_edges(2:end));
    lat_bin = 0.5*(lat_edges(1:end-1) + lat_edges(2:end));
    
    figure(3);
    subplot(1,2,1)
    tmp=histogram2(LONs(selmask2), LATs(selmask2),lon_edges, lat_edges,'DisplayStyle','tile','ShowEmptyBins','off','normalization','pdf');
    wanom_pdf{tt} = tmp.Values';
    colorbar;
    
    subplot(1,2,2)
    tmp=histogram2(LONs(selmask), LATs(selmask),lon_edges, lat_edges,'DisplayStyle','tile','ShowEmptyBins','off','normalization','pdf');
    canom_pdf{tt} = tmp.Values';
    colorbar
    
    maxpdf_val = max([max(wanom_pdf{tt}(:)), max(canom_pdf{tt}(:))]);
    pdf_levs{tt} = [0.1:0.2:1]*maxpdf_val;
    
    
    %% save the following variables for plotting:
    % mean SST:
    % selected locations of warm and cold fronts.
    
    save('locationsPDF_of_extreme_effective_SSTgrad_more_time_consistent.mat','SST_LS', 'wanom_pdf', 'canom_pdf','pdf_levs','lon_bin','lat_bin');
    
end

cmap = getPyPlot_cMap('coolwarm',64);
pyRdBu = getPyPlot_cMap('RdBu_r');
redblue_all = pyRdBu ;
redblue_mid = pyRdBu(33:2:96,:);

figure
pcolor(LON_sub, LAT_sub,mean(CFRC.(CN),3,'omitnan')); shading flat;
%colorbar; caxis([-4, 4]*10^(-4))
if strcmp(CN, 'flowers')
    colormap(redblue_all);
    caxis([-1.5,1.5]);
else
    colormap(redblue_mid);
    caxis([-0.5, 0.5]);
end
hold on
hb =colorbar; %caxis([-2.5, 2.5]*10^(-4))
set(get(hb,'xlabel'), 'string', '\bf{u} \cdot \nablaSST (K/hr)');
[C,h]= contour(LON_sub, LAT_sub, mean(SST_LS.(CN),3,'omitnan')-273.15, [25:0.5:28],'--k');
clabel(C,h, [25:0.5:28]);
% lon_bin = -59+0.25/2:0.25:-48-0.25/2;
% lat_bin = 8+0.25/2: 0.25: 18-0.25/2;
contour(lon_bin, lat_bin,canom_pdf,pdf_levs, 'b');
contour(lon_bin, lat_bin,wanom_pdf,pdf_levs, 'm');

% plot(tail_aTend_locs.lon, tail_aTend_locs.lat,'.c', 'linewidth',0.2);
% plot(LONs(selmask2), LATs(selmask2), '.m', 'linewidth',0.2);
% contour(LONs, LATs, std(aTtend.(CN),1,3),[1:0.5:2.5].*10^(-4),'k');

%plot(LONs(selmask2), LAT(selmask2),'.c');
%plot(LONs(selmask3), LAT(selmask3),'.y');
title([CN '-favored day averaged'])

xlabel('longitude'); ylabel('latitude');
set(gca,'fontsize',14);
axis('equal')

%%%%%%%%% obsolete: %%%%%%%%%%%%%%
for i=1:1
    
    % This is an interesting plot:
    aTbin = 0.5*(aTedges(1:end-1) + aTedges(2:end));
    lncols = getPyPlot_cMap('Paired',10);
    figure(7); clf;
    for i = 1:4
        
        subplot(2,2,i)
        num_of_occur_per_cloudinessbin = sum(occur_counts{i},1);
        %bids = find(num_of_occur_per_cloudinessbin>0);
        num_of_occur_per_cloudinessbin(num_of_occur_per_cloudinessbin==0)=NaN;
        %normalized_occur = [];
        % for ib =1: length(bids)
        normalized_occur = occur_counts{i}./num_of_occur_per_cloudinessbin;
        %end
        clear lgdstr h
        cnt = 0;
        for il = 1:size(normalized_occur,2)
            if ~isnan(num_of_occur_per_cloudinessbin(il))
                cnt = cnt+1;
                lgdstr{cnt} = [num2str(CFedges(il)) '~' num2str(CFedges(il+1))];
                h(cnt) = plot(aTbin,  normalized_occur(:,il), 'color', lncols(il,:), 'linewidth', 1.2);
                hold on;
                set(gca,'yscale','log');
                grid on
            end
        end
        plot([0 0], [10^(-5), 10^0],'--k');
        hlgd=legend(h, lgdstr);
        set(hlgd,'fontsize',11,'Location','northwest','Orientation','vertical');
        xlim([-4 4]*10^(-4));
        ylim(10.^[-5,0]);
        xlabel('\bf{u} \cdot \nablaSST (K/hr)')
        ylabel({'normalized occurence', 'for each cloudiness bin'})
        set(gca,'fontsize',13);
        title([cloudtypes{i} '-favored regime'],'fontsize',13);
        
    end
    figname = 'probability_of_airTtend_at_different_cloudiness_level.jpg';
    xc_savefig(gcf, './Figs/probability_exploration', figname, [0 0 10 8]);
    
    
    
    % do this analysis for laplacian:
    lncols = getPyPlot_cMap('Paired',10);
    SSTLpbin = 0.5*(SSTLpEdges(1:end-1) + SSTLpEdges(2:end));
    figure(8); clf;
    for i = 1:4
        subplot(2,2,i)
        num_of_occur_per_cloudinessbin = sum(occur_counts_laplacian{i},1);
        %bids = find(num_of_occur_per_cloudinessbin>0);
        num_of_occur_per_cloudinessbin(num_of_occur_per_cloudinessbin==0)=NaN;
        %normalized_occur = [];
        % for ib =1: length(bids)
        normalized_occur = occur_counts_laplacian{i}./num_of_occur_per_cloudinessbin;
        %end
        clear lgdstr h
        cnt = 0;
        for il = 1:size(normalized_occur,2)
            if ~isnan(num_of_occur_per_cloudinessbin(il))
                cnt = cnt+1;
                lgdstr{cnt} = [num2str(CFedges(il)) '~' num2str(CFedges(il+1))];
                h(cnt) = plot(SSTLpbin*10^2,  normalized_occur(:,il), 'color', lncols(il,:), 'linewidth', 1.2);
                hold on;
                set(gca,'yscale','log');
                grid on
            end
        end
        plot([0 0], [10^(-5), 10^0],'--k');
        hlgd=legend(h, lgdstr);
        set(hlgd,'fontsize',11,'Location','northeastout','Orientation','vertical');
        xlim([-1.6 1.6]*10^(-9+2));
        ylim(10.^[-5,0]);
        xlabel('{\nabla}^2SST (K/(10km)^2)')
        ylabel({'normalized occurence', 'for each cloudiness bin'})
        set(gca,'fontsize',13);
        title([cloudtypes{i} '-favored regime'],'fontsize',13);
        
    end
    figname = 'probability_of_SSTLaplacian_at_different_cloudiness_level.jpg';
    xc_savefig(gcf, './Figs/probability_exploration', figname, [0 0 12 8]);
    
end