% This script is used to compare the satellite SST and the SST sampled by
% Ron Brown, on days that have moderate SST variability. 

% satellite SST is interpolated linearly in space and time to form
% comparison with the RHB SST.

% according to Gary,  so days are better than others. (So I will need to be
% careful about it. (why though?, is it because of the spatial scale??)
%

% load in the RHB data:
RHB_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
load([RHB_datadir filesep 'rhb_daily_grouped_10min_data_0909latest.mat']);
% load in the coastline I had:
addpath('/Users/xchen/Documents/MATLAB/plotting_resources/gshhs_coastline_data');
load(['NACoast_shoreline_h.mat']);   % In this version of the coast line, the Barbados is missing..

% RHB_ds:
time_mask =RHB_ds.time < datenum(2020,1,10);
RHB_lon = RHB_ds.lon(time_mask);
RHB_lat = RHB_ds.lat(time_mask);

strong_SSTvar_dates = rhbdates.strong_SSTvar;
moderate_SSTvar_dates = rhbdates.moderate_SSTvar;

all_DOYs = 1+[strong_SSTvar_dates moderate_SSTvar_dates]-basetimenum;
all_DOYs = sort(all_DOYs);

for id = 1: length(DOYs_extra) %length(all_DOYs)
    % DOY = all_DOYs(id);
    DOY = DOYs_extra(id);
    %path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];
    %dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');
    dateID_mmmdd = datestr(basetimenum + (DOY-1),'mmmdd');
    load(['GOES_SST_' dateID_mmmdd '.mat']);
    disp(['working on ' dataID_date]);

  % load in data on a specific day:

    % make plots to check the extraction:
    QC_mask = ones(size(GOES_ATOMIC.quality_level));
    QC_mask(GOES_ATOMIC.quality_level~=5) = NaN;
    masked_SST = GOES_ATOMIC.sea_surface_temperature.*QC_mask;
    
    
    % RHB_ds:
    time_mask =(RHB_ds.time>basetimenum+(DOY-1)) & (RHB_ds.time<basetimenum+DOY);
    RHB_lon = RHB_ds.lon(time_mask);
    RHB_lat = RHB_ds.lat(time_mask);
    
    RHB_time = RHB_ds.time(time_mask);
    
    hfig(1) = figure(1); clf;
    for it =1:nhr
        %subplot(2,1,1)
        pcolor(GOES_ATOMIC.lon(:,it), GOES_ATOMIC.lat(:,it), squeeze(GOES_ATOMIC.sea_surface_temperature(:,:,it))'-273.15);
        shading interp;
        hold on;
        hb = colorbar;
        caxis([25 28.5])
        % add coastline:
        plot(coast.lon, coast.lat, '-k', 'linewidth',2);
        % add RHB track on this selected day:
        plot(RHB_lon,RHB_lat,'-r','linewidth',1.5);
        title(datestr(GOES_ATOMIC.time_num(it),'mmm-dd-yyyy HH:MM'),'(GOES-16 2km SST)')
        set(gca,'fontsize',15,'TickDir','both');
        hold off
        pause(0.5);
        frame = getframe(hfig(1));
        im{it} = frame2im(frame);
    end
   
    VideoName=[dateID_mmmdd '_GOES16_SST_L3C_2km_hrly.gif'];
    
    for it = 1:length(im)
        [A,map] = rgb2ind(im{it},256);
        if it == 1
            imwrite(A,map,VideoName,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,VideoName,'gif','WriteMode','append','DelayTime',1);
        end
    end
    
    % so now, we can do some interpolation and compare the SST:
  lon_sat = double(GOES_ATOMIC.lon(:,1));
  lat_sat = double(GOES_ATOMIC.lat(:,1));
  time_sat = GOES_ATOMIC.time_num;
  SST_sat = permute(GOES_ATOMIC.sea_surface_temperature,[2,1,3]);
  
  
  
  
  satSST_atRHB = interp3(lon_sat, lat_sat,time_sat, SST_sat, RHB_lon, RHB_lat, RHB_time);

    
    
  % quickly check results:
  figure(10); clf;
  offset=-273.15;
  plot(RHB_time, satSST_atRHB+offset,'.-k');  
  hold on
  plot(RHB_time, RHB_ds.tsea(time_mask), '.-r');
  %plot(RHB_time, RHB_ds.tskin(time_mask), '.-b');
  datetick('x')
  hold off
  xlabel('time')
  ylabel('SST (Celcius)');
  hlgd = legend({'satellite','RHB tsea', 'RHB tskin'});
  set(hlgd,'location','northwest');
  set(gca,'fontsize',14);
  title(dateID_mmmdd);
  xc_savefig(gcf,'./',[dateID_mmmdd '_GOES16_versus_RHB_tsea.jpg'],[0 0 8 6]);
  
  
  G16_SST_collocated.(dateID_mmmdd) = satSST_atRHB;
end

% compare the SST variability from two sources;


