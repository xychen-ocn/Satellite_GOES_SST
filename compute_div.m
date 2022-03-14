function div = compute_div(lon, lat, u, v)
% du/dx + dv/dy
[LON, LAT] =meshgrid(lon, lat);

latmid = 0.5*(lat(1:end-1)+lat(2:end));
lonmid = 0.5*(lon(1:end-1)+lon(2:end));
lon2x = 111E3.*cosd(LAT);
lat2x = 111E3;

% forward difference:
du = diff(u,1,2);
dx = diff(LON,1,2).* lon2x(:,1:end-1);

dv = diff(v,1,1);
dy = diff(LAT,1,1).* lat2x;

% interpolate back to the original grid:
[LON1, LAT1] = meshgrid(lonmid, lat);
du_dx = interp2(LON1, LAT1, du./dx, LON, LAT);

[LON2, LAT2] = meshgrid(lon, latmid);
dv_dy = interp2(LON2, LAT2, dv./dy, LON, LAT);

div = du_dx + dv_dy;

return