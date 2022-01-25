% scripts used for satellite data:
clear all; clc; close all;

% 1. read in the L3C raw data
% set up data directory.
basetimenum = datenum('2020-01-01','yyyy-mm-dd');
data_root = '/Volumes/g16';
DOY = 9;
path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];

dataFN_suffix = '-STAR-L3C_GHRSST-SSTsubskin-ABI_G16-ACSPO_V2.70-v02.0-fv01.0.nc';
dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');

hourvec = 0:23;
nhr = length(hourvec);

ATOMIC_area.lon = [-62 -48];
ATOMIC_area.lat = [8 18];

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

for id = 1:length(all_DOYs)
    DOY = all_DOYs(id);
    path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];
    dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');
    disp(['working on ' dataID_date]);
    for it = 1 :nhr
        % dynamically construct data file:
        dataID_time = [num2str(hourvec(it),'%2.2i') '0000'];
        dataFN = [dataID_date dataID_time dataFN_suffix];
        
        % read in the data:
        absFN = [path2data filesep dataFN];
        
        % if extra argument present, then this matlab function will take data
        % from a subset of the netCDF file.
        
        GOES_SST =read_netCDF_into_matlab_structure(absFN);
        if mod(it,6)
            disp(['finished reading in ' num2str(it) 'hrs..']);
        end
        if it ==1
            [tmp, lon_stid]=min(abs(GOES_SST.lon-ATOMIC_area.lon(1)));
            [tmp, lon_edid]=min(abs(GOES_SST.lon-ATOMIC_area.lon(2)));
            
            [tmp, lat_stid]=min(abs(GOES_SST.lat-ATOMIC_area.lat(1)));
            [tmp, lat_edid]=min(abs(GOES_SST.lat-ATOMIC_area.lat(2)));
            if lat_stid < lat_edid
                lat_inc = 1;
            else
                lat_inc = -1;
            end
            
            
            chunk.Nlon = lon_edid - lon_stid +1;
            chunk.Nlat = lat_edid - lat_stid +1;
        end
        % taking subset of the data manually:
        datafields = fieldnames(GOES_SST);
        GOES_ATOMIC.lon(:,it) = GOES_SST.lon(lon_stid:lon_edid);
        GOES_ATOMIC.lat(:,it) = GOES_SST.lat(lat_stid:lat_inc:lat_edid);
        GOES_ATOMIC.time(it) = GOES_SST.time;
        GOES_ATOMIC.time_num(it) = GOES_SST.time_num;
        GOES_ATOMIC.crs(it) = GOES_SST.crs;
        for j = 1:length(datafields)
            FN = datafields{j};
            if j>3 && j<14
                
                GOES_ATOMIC.(FN)(:,:,it) = GOES_SST.(FN)(lon_stid:lon_edid, lat_stid:lat_inc:lat_edid);
                
            end
        end
        
    end
    save(['GOES_SST_' datestr(basetimenum + (DOY-1),'mmmdd') '.mat'],'GOES_ATOMIC');
end
%     else
%     
%     % take chunk of data out of the large netCDF
%     chunk.lon_stid = lon_stid;
%     chunk.lat_stid = min([lat_stid, lat_edid]);
%     GOES_ATOMIC_tmp = read_netCDF_into_matlab_structure(absFN, chunk);
%     
    % find the SST along the RHB trajectory. (spatial, temporal
    % interpolation (linear?)
    
    
    % make plots to check the extraction:
    QC_mask = ones(size(GOES_ATOMIC.quality_level));
    QC_mask(GOES_ATOMIC.quality_level~=5) = NaN;
    masked_SST = GOES_ATOMIC.sea_surface_temperature.*QC_mask;
    
    
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
   
    VideoName='Jan09_GOES16_SST_L3C_2km_hrly.gif';
    
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
  RHB_time = RHB_ds.time(time_mask);
 
  methods = {'linear','cubic','spline'};
  for i = 1:length(methods)
      MN = methods{i};
      satSST_atRHB = interp3(lon_sat, lat_sat,time_sat, SST_sat, RHB_lon, RHB_lat, RHB_time);
  end
    
    
  % quickly check results:
  figure(10);
  offset=-273.15;
  plot(RHB_time, satSST_atRHB+offset,'.-k');  
  hold on
  plot(RHB_time, RHB_ds.tsea(time_mask), '.-r');
  plot(RHB_time, RHB_ds.tskin(time_mask), '.-b');
  datetick('x')
  hold off
  xlabel('time')
  ylabel('SST (Celcius)');
  legend({'satellite','RHB tsea', 'RHB tskin'})
  set(gca,'fontsize',14);
  title('Jan-09-2020');
  xc_savefig(gcf,'./','Jan09-quick_check_promising.jpg',[0 0 8 6]);
  





% write a function to read in the parameter of interest from the raw data;

% 2. compute SST laplacian at each pixel (follow the method in Li and
% Carbone 2012 method.)


% 3. for each pixel, record the SST laplacian before cloud masks the pixel
% and after the pixel first clear up. 

% 4. then check how the distribution looks like.