function download_GOES_SSTdata(days_of_year, datasvdir, varargin)
%

switch nargin
    case(2)
        requested_area.lon = [-62 -48];
        requested_area.lat = [8 18];
        
    case(3)
        requested_area = varargin;
end
% 1. read in the L3C raw data
% set up data directory.
basetimenum = datenum('2020-01-01','yyyy-mm-dd');
data_root = '/Volumes/g16';

dataFN_suffix = '-STAR-L3C_GHRSST-SSTsubskin-ABI_G16-ACSPO_V2.70-v02.0-fv01.0.nc';

hourvec = 0:23;
nhr = length(hourvec);


for id = 1:length(days_of_year)
   
    DOY = days_of_year(id);
    path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];
    dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');
    
    
    matDataFN = ['GOES_SST_' datestr(basetimenum + (DOY-1),'mmmdd') '.mat'];
    
    % download if local data copy not exist:
    if exist([datasvdir filesep matDataFN], 'file')==0
        
        disp(['downloading ' dataID_date]);
        
   % for it = 1 :nhr
   ncfiles = dir([path2data filesep '*.nc']);
   dataFNs = {ncfiles.name};
   
   for it = 1:length(dataFNs)
        % dynamically construct data file:
        %dataID_time = [num2str(hourvec(it),'%2.2i') '0000'];
        %dataFN = [dataID_date dataID_time dataFN_suffix];
        dataFN = dataFNs{it};
        
        % read in the data:
        absFN = [path2data filesep dataFN];
        
        % if extra argument present, then this matlab function will take data
        % from a subset of the netCDF file.
        
        GOES_SST =read_netCDF_into_matlab_structure(absFN);
        if mod(it,6)==0
            disp(['finished reading in ' num2str(it) 'hrs..']);
        end
        if it ==1
            [tmp, lon_stid]=min(abs(GOES_SST.lon-requested_area.lon(1)));
            [tmp, lon_edid]=min(abs(GOES_SST.lon-requested_area.lon(2)));
            
            [tmp, lat_stid]=min(abs(GOES_SST.lat-requested_area.lat(1)));
            [tmp, lat_edid]=min(abs(GOES_SST.lat-requested_area.lat(2)));
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
                % can be potentially improved here.
                GOES_ATOMIC.(FN)(:,:,it) = GOES_SST.(FN)(lon_stid:lon_edid, lat_stid:lat_inc:lat_edid);
                
            end
        end
        
    end
    % save data:
    save([datasvdir filesep matDataFN],'GOES_ATOMIC');
    
    end
end

return