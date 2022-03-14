function sa = get_satellite_SST_anom_along_RHBtraj(SSTa_maps, LON, LAT, times, RHB_full)
% the input RHB_full contains UTC time;
% interpolate the SSTa data onto the RHB sample time (every 10-min) for the
% day of interest. the day is the day of the input "times".
%
tday = unique(floor(times));

time_mask = RHB_full.time>tday &  RHB_full.time<tday+1;
RHB_lon = RHB_full.lon(time_mask);
RHB_lat = RHB_full.lat(time_mask);
RHB_time = RHB_full.time(time_mask);

sa.values = interp3(LON(1,:), LAT(:,1), times, SSTa_maps, RHB_lon, RHB_lat, RHB_time);
sa.time = RHB_time;
sa.lon = RHB_lon;
sa.lat = RHB_lat;



return