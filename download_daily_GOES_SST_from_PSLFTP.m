% This script is used to download satellite data from the ftp;
% (I already downloaded all the days I need, this is script below is saved
% in case in need of it further.)

% scripts used for satellite data:
clear all; clc; close all;

% 1. read in the L3C raw data
% set up data directory.
basetimenum = datenum('2020-01-01','yyyy-mm-dd');
data_root = '/Volumes/g16';

dataFN_suffix = '-STAR-L3C_GHRSST-SSTsubskin-ABI_G16-ACSPO_V2.70-v02.0-fv01.0.nc';

hourvec = 0:23;
nhr = length(hourvec);

ATOMIC_area.lon = [-62 -48];
ATOMIC_area.lat = [8 18];

% load in the RHB data:
RHB_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/data';
load([RHB_datadir filesep 'rhb_daily_grouped_10min_data_0909latest.mat'],'rhbdates');


strong_SSTvar_dates = rhbdates.strong_SSTvar;
moderate_SSTvar_dates = rhbdates.moderate_SSTvar;

all_DOYs = 1+[strong_SSTvar_dates moderate_SSTvar_dates]-basetimenum;
all_DOYs = sort(all_DOYs);

DOYs_extra =1+[datenum(2020,1,19), datenum(2020, 2,7), datenum(2020, 2, 12)] - basetimenum;

for id = 1:length(DOYs_extra) %length(all_DOYs)
   % DOY = all_DOYs(id);
    DOY = DOYs_extra(id);
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