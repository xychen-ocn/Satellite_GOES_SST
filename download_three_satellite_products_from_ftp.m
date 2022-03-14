% This script will download subsets of satellite SST data to local
% directory, for three different type of datasets.

clear all; clc;

local_root = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
DatasetNames ={'g16','gpblend','ostia'};
numDS = length(DatasetNames);

% set up storage place locally;
for i = 1:numDS
    DN = DatasetNames{i};
    datadir_tmp = [local_root filesep DN];
    if ~exist(datadir_tmp)
        mkdir(datadir_tmp);
    end
    local_datadir.(DN) = datadir_tmp;
end

% specify region of interests:
atomic_area.lon = [-62, -42];
atomic_area.lat = [5, 25];

% start downloading data using functions:
% 1. for g16:

%days_of_year = 1:59;
days_of_year = 334:365;
download_GOES_SSTdata(days_of_year, local_datadir.g16, 'subset_area',atomic_area,'dataroot','/Volumes/2019');   % somehow the ftp data can't be open?

% 2. for gpblend and ostia:
for i = 2:numDS
    DN = DatasetNames{i};
    download_L4_products(DN, local_datadir.(DN), 'subset_area',atomic_area); % 'saveMode','Together');
end



