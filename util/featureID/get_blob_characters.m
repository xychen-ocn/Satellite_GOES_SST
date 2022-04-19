function blob_stats = get_blob_characters(blob_samples, spatialres)
% purpose: this script is aimed to get out the following information from
% the structure above and save it into a better format
% Input: blob_samples (matlab structure of length nblob)
%     
% Output: blob_stats (structure 1x1)
%   with field:  SST_anom (strength), array of nblob x1
%             :  SST_background
%             :  EqvDiam  (size)
%             :  Area (size)
%             :  MajAxisLen, MinAxisLen (shape), array of nblob x 2
%             :  Maj/Min Ratio (shape) 
%             :  Orientation (in case needed)
% Note: somehow, maxFeretCoord has been recorded in...
%nblob = length(blob_samples);

blob_stats.SST_anom = [blob_samples.max_SSTa];
blob_stats.SST_bg = [blob_samples.ave_SSTbg];
stats = [blob_samples.stats_selected];
blob_stats.EqvDiam_km = [stats.EqvDiam].* spatialres;
blob_stats.Area_kmsq = [stats.Area].* spatialres.^2;
blob_stats.MajAxisLen = [stats.MajAxisLen];
blob_stats.MinAxisLen = [stats.MinAxisLen];
blob_stats.AxisRatio = blob_stats.MajAxisLen./blob_stats.MinAxisLen;
blob_stats.Orientation = [stats.orientation];


return