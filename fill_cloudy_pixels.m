function [SST_map_filled, filled_flag] = fill_cloudy_pixels(SSTL3, LON, LAT, cld_thres)
% use linear interpolation to fill data gap:
% input: SSTL3 (Nlat, Nlon, Nt)
%  LON (Nlat, Nlon), LAT 

nt = size(SSTL3,3);
disp(['ntime=' num2str(nt)]);
SST_map_filled = zeros(size(SSTL3));
filled_flag = false(nt,1);

for it = 1:nt
    SST_map = squeeze(SSTL3(:,:,it));
    SST_map_copy = SST_map;
    cloud_flag = isnan(SST_map);
    
    total_pixelnum = numel(SST_map);
    cloudy_pixelnum = length(find(cloud_flag));
    cloud_pcen = cloudy_pixelnum/total_pixelnum;
    
    mean_SST = mean(SST_map,[1,2],'omitnan');
    if ~isnan(mean_SST) && cloud_pcen<=cld_thres 
        % not every pixel is cloud and cloud pcen< a certain threshold
        % (filled the data)
        % get the interpolation operator from non-cloudy pixels:
        SST_intpFunc = scatteredInterpolant(LON(~cloud_flag), LAT(~cloud_flag), SST_map(~cloud_flag),'natural', 'nearest');
        SST_map_copy(cloud_flag)=SST_intpFunc(LON(cloud_flag), LAT(cloud_flag));
        
        SST_map_filled(:,:,it) = SST_map_copy;
        filled_flag(it) = true;
    else
        % do not fill the SST map.
        SST_map_filled(:,:,it) = SST_map;
    end
       
    if mod(round(it/24*100),10)==0
        disp(['completed ' ,num2str(round(it/24*100)) '% of inputs..']);
    end
end

return