function download_L4_products(dsName, datasvdir, varargin)
% this function is used to download two L4 products:
%
default_area.lon = [-62 -48];
default_area.lat = [8 18];

default_root = ['/Volumes/sst_raw'];


p = inputParser;
addRequired(p, 'dsName', @ischar);
addRequired(p, 'datasvdir', @ischar);
addParameter(p, 'subset_area', default_area, @isstruct);
addParameter(p, 'dataroot', default_root, @ischar);


parse(p, dsName, datasvdir,varargin{:});
dsName = p.Results.dsName;
datasvdir = p.Results.datasvdir;
requested_area = p.Results.subset_area;
data_root = p.Results.dataroot;

path2data = [data_root filesep dsName];

if strcmpi(dsName, 'gpblend')
    dataFN_suffix = '-OSPO-L4_GHRSST-SSTfnd-Geo_Polar_Blended-GLOB-v02.0-fv01.0.nc';
    
elseif strcmpi(dsName, 'ostia')
    dataFN_suffix = '-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc';
    
end


files = dir([path2data filesep '*' dataFN_suffix]);
filenames = {files.name};

nfiles = length(filenames);
for i = 1:nfiles
    absFN = [path2data filesep filenames{i}];
    
    if i ==1
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
        
        SSTdata =read_netCDF_into_matlab_structure(absFN, 'starts', starts, 'counts', counts);
        
        if i ==1
            disp(['check lon0:' num2str(SSTdata.lon(1)) ';\n lonN:' num2str(SSTdata.lon(end))]);
            disp(['check lat0:' num2str(SSTdata.lat(1)) ';\n latN:' num2str(SSTdata.lat(end))]);
        end

        
        data.lon(:,i) = SSTdata.lon;
        if lat_inc == -1
            SSTdata.lat = flipud(SSTdata.lat);
        end
        data.lat(:,i) = SSTdata.lat;
        data.time(i) = SSTdata.time;
        data.time_num(i) = SSTdata.time_num;
        
        datafields = fieldnames(SSTdata);
        % need to double check here:
        for j = 1:length(datafields)
            FN = datafields{j};
            if length(size(SSTdata.(FN)))==2 && numel(SSTdata.(FN))>length(SSTdata.(FN))
                if lat_inc == -1
                    data.(FN)(:,:,i) = SSTdata.(FN)(:, end:-1:1);
                else
                    data.(FN)(:,:,i) = SSTdata.(FN);
                end
                
            end
        end
        
     if mod(round(i/nfiles*100) , 10)==0
         disp(['finished reading ' num2str(round(i/nfiles*100)) '%...']);
     end
end

% eval([upper(dsName) '_ATOMIC = data;']);
% clear data
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


matDataFN = [dsName '_L4_daily_JanFeb_ATOMIC_lon' lonstr{1}, '_TO_' lonstr{2}  ...
                '_lat' latstr{1}, '_TO_' latstr{2} '.mat'];
save([datasvdir filesep matDataFN], 'data');






end