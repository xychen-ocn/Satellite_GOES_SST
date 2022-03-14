function output = compute_cloudfreq_and_anom(SSTmaps, mean_CF)
% compute cloudiness freq from input SST maps and compute anomalies
% relative to the input mean state.

CF_flag = isnan(SSTmaps);

CF = sum(double(CF_flag),3)./size(SSTmaps,3);
CF_anom= (CF - mean_CF)./mean_CF;

output.CF = CF;
output.CF_anom = CF_anom;
output.SSTmaps = SSTmaps;


end