function [SST_map_filled, filled_flag] = fill_cloudy_pixels(SSTL3, LON, LAT, time, cld_thres)
% use linear interpolation to fill data gap:
% input: SSTL3 (Nlat, Nlon, Nt)
%  LON (Nlat, Nlon), LAT 
% date: updated on Jan 26, 2022 (interpolation in time at each grid point
% with intermittent cloud contamination (20% of the time) first, and then do spatial interpolation)
%

[NY, NX, nt] = size(SSTL3);
%nt = size(SSTL3,3);
disp(['ntime=' num2str(nt)]);
SST_map_filled = zeros(size(SSTL3));
filled_flag = false(nt,1);

% use temporal interpolation first to fill out some intermittent gaps:
% then use spatial inteprolation to fill out the rest..
thr = (time-time(1))*24;
        
for j = 1:NY
    for i = 1:NX
        SST_ts = squeeze(SSTL3(j,i,:));
        
        % find valid data:
        valid = ~isnan(SST_ts);
        good_data_prc = round(length(find(valid==1))/length(SST_ts)*10)/10; % good data percentage
        
        if good_data_prc>=0.6
            SST_ts_intp(j,i,:) = interp1(thr(valid), SST_ts(valid), thr,'linear');
%             plot(thr, SST_ts_intp,'.-');
%             pause
        
        else
            SST_ts_intp(j,i,:) = SST_ts;    % no filling.
        end
        
        %pause
    end
end

% spatial intepolration:
for it = 1:nt
    SST_map = squeeze(SST_ts_intp(:,:,it));
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