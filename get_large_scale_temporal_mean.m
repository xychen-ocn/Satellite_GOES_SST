function SST_mean = get_large_scale_temporal_mean(dsName, wz, LPF_scale_km, region)
% Purpose: get the saved large scale temporal mean data requested.
% Inputs: dsName (char)
%         wz: window size (of averaging)
%         LPF_scale_km: (cutoff scale in km for low-pass filtering)
%         regino (2 x 2 array) to get subset of data out ([lon0, lonN;
%         lat0, latN])
% Outputs: 
%         SST_LS (structure with field: values and time, aveMethod);
%
% Date: Jan 30, 2022 (XYC)
%

dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';

datadir = [dataroot filesep lower(dsName)];       % in case the dsName is not in lowercase;

matFN = ['TemporalAveraged_LowPassFiltered_thres' num2str(LPF_scale_km) 'km_SST.mat'];

fieldN = ['wndsz_' num2str(wz) 'd'];


if wz < 15
    % moving average:
    load([datadir filesep matFN], [upper(dsName) '_movaved']);
    %feval('load',[datadir filesep matFN], [upper(dsName) '_movaved']);
    eval(['data = ' upper(dsName) '_movaved;']);
    SST_mean.aveMethod = 'moving';
    
else
    % fixed average:
    load([datadir filesep matFN], [upper(dsName) '_aved']);
    eval(['data = ' upper(dsName) '_aved;']);
    SST_mean.aveMethod = 'fixed';
    
end

% get subset mask:
lonmask = data.lon>=region(1,1) & data.lon<=region(1,2);
latmask = data.lat>=region(2,1) & data.lat<=region(2,2);

SST_mean.values = permute(data.(fieldN).SST(lonmask, latmask,:), [2,1,3]);
SST_mean.LPF = data.(fieldN).SST_LS(latmask, lonmask,:);
SST_mean.time = data.(fieldN).time;
SST_mean.lon = data.lon(lonmask);
SST_mean.lat = data.lat(latmask);



return