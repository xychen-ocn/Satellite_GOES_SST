function download_GOES_SSTdata(days_of_year, datasvdir, varargin)
%
default_area.lon = [-62 -48];
default_area.lat = [8 18];

default_root = '/Volumes/sst_raw/g16';


p = inputParser;
addRequired(p, 'days_of_year', @isnumeric);
addRequired(p, 'datasvdir', @ischar);
addParameter(p, 'subset_area', default_area, @isstruct);
addParameter(p, 'dataroot', default_root, @ischar);

parse(p, days_of_year, datasvdir,varargin{:});
days_of_year = p.Results.days_of_year;
datasvdir = p.Results.datasvdir;
requested_area = p.Results.subset_area;
data_root = p.Results.dataroot;


% switch nargin
%     case(2)
%         requested_area.lon = [-62 -48];
%         requested_area.lat = [8 18];
%         
%     case(3)
%         requested_area = varargin;
% end
% 1. read in the L3C raw data
% set up data directory.
%basetimenum = datenum('2020-01-01','yyyy-mm-dd');
basetimenum = datenum('2019-01-01','yyyy-mm-dd');


dataFN_suffix = '-STAR-L3C_GHRSST-SSTsubskin-ABI_G16-ACSPO_V2.70-v02.0-fv01.0.nc';

%hourvec = 0:23;
%nhr = length(hourvec);

if requested_area.lon(1)<0
    lonstr{1} = [num2str(-requested_area.lon(1)), 'W'];
    lonstr{2} = [num2str(-requested_area.lon(2)), 'W'];
else
    lonstr{1} = [num2str(requested_area.lon(1)), 'E'];
    lonstr{2} = [num2str(requested_area.lon(2)), 'E'];
end

if requested_area.lat(1)>0
    latstr{1} = [num2str(requested_area.lat(1)), 'N'];
    latstr{2} = [num2str(requested_area.lat(2)), 'N'];
else
    latstr{1} = [num2str(-requested_area.lat(1)), 'S'];
    latstr{2} = [num2str(-requested_area.lat(2)), 'S'];
end



for id = 1:length(days_of_year)
   
    DOY = days_of_year(id);
    path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];
    dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');
    
    
    matDataFN = ['GOES_SST_' num2str(DOY,'%3.3d') '_' datestr(basetimenum + (DOY-1),'mmmdd') ...
                '_lon' lonstr{1}, '_TO_' lonstr{2}  ...
                '_lat' latstr{1}, '_TO_' latstr{2} '.mat'];
    
    % download if local data copy not exist:
    if exist([datasvdir filesep matDataFN], 'file')==0
        
        disp(['downloading ' dataID_date]);
        
   % for it = 1 :nhr
   ncfiles = dir([path2data filesep '*.nc']);
   dataFNs = {ncfiles.name};
   
   clear GOES_ATOMIC
   disp(['nt = ' num2str(length(dataFNs))]);
   for it = 1:length(dataFNs)
        % dynamically construct data file:
        %dataID_time = [num2str(hourvec(it),'%2.2i') '0000'];
        %dataFN = [dataID_date dataID_time dataFN_suffix];
        dataFN = dataFNs{it};
        
        % read in the data:
        absFN = [path2data filesep dataFN];
        
        % if extra argument present, then this matlab function will take data
        % from a subset of the netCDF file.
        
        % read lon and lat out to get indices to get subset region out:
        
        if it ==1
            lon = ncread(absFN, 'lon');
            lat = ncread(absFN, 'lat');
            
            [tmp, lon_stid]=min(abs(lon-requested_area.lon(1)));
            [tmp, lon_edid]=min(abs(lon-requested_area.lon(2)));
            
            [tmp, lat_stid]=min(abs(lat-requested_area.lat(1)));
            [tmp, lat_edid]=min(abs(lat-requested_area.lat(2)));
            
            chunk.Nlon = lon_edid - lon_stid +1;
            
            
            if lat_stid < lat_edid
                lat_inc = 1;
                chunk.Nlat = lat_edid - lat_stid +1;
                starts = [lon_stid, lat_stid, 1];
            else
                lat_inc = -1;
                chunk.Nlat = lat_stid - lat_edid +1;
                starts = [lon_stid, lat_edid, 1];
            end
            
            
        end
        
        
        counts = [chunk.Nlon, chunk.Nlat, Inf];
        
        GOES_SST =read_netCDF_into_matlab_structure(absFN, 'starts', starts, 'counts', counts);
        
        if it ==1
            disp(['check lon0:' num2str(GOES_SST.lon(1)) ';\n lonN:' num2str(GOES_SST.lon(end))]);
            disp(['check lat0:' num2str(GOES_SST.lat(1)) ';\n latN:' num2str(GOES_SST.lat(end))]);
        end
        if mod(it,6)==0
            disp(['finished reading in ' num2str(it) 'hrs..']);
        end

        % taking subset of the data manually:
        datafields = fieldnames(GOES_SST);
        GOES_ATOMIC.lon(:,it) = GOES_SST.lon;
        %[ascend_lat, sid] = sort(GOES_SST.lat,'ascend');\
        if lat_inc == -1
            GOES_SST.lat = flipud(GOES_SST.lat);
        end
        GOES_ATOMIC.lat(:,it) = GOES_SST.lat;
        GOES_ATOMIC.time(it) = GOES_SST.time;
        GOES_ATOMIC.time_num(it) = GOES_SST.time_num;
        %GOES_ATOMIC.crs(it) = GOES_SST.crs;
        for j = 1:length(datafields)
            FN = datafields{j};
            if length(size(GOES_SST.(FN)))==2 && numel(GOES_SST.(FN))> length(GOES_SST.(FN))
                % can be potentially improved here.
                if lat_inc == -1
                    GOES_ATOMIC.(FN)(:,:,it) = GOES_SST.(FN)(:, end:-1:1);
                else
                    GOES_ATOMIC.(FN)(:,:,it) = GOES_SST.(FN);
                end
                
            end
        end
        
    end
    % save data:
    save([datasvdir filesep matDataFN],'GOES_ATOMIC');
    
    end
end

return