function [day, night]=group_data_into_daytime_and_nighttime_samples(SSTmaps, time)
% separate input maps into day and nights:

time_dailyhrs = (time - floor(time))*24;
nightTime_mask = false(size(time_dailyhrs));
nightTime_mask(time_dailyhrs<=10.0)=true;
nightTime_mask(time_dailyhrs>21.0)=true;

dims = size(SSTmaps);
n = length(dims);
if n==3
    night.SSTmaps = SSTmaps(:,:,nightTime_mask);
    day.SSTmaps = SSTmaps(:,:,~nightTime_mask);
    
elseif n==2
    night.SSTmaps = SSTmaps(:, nightTime_mask);
    day.SSTmaps = SSTmaps(:, ~nightTime_mask);
end

night.time = time(nightTime_mask);
day.time = time(~nightTime_mask);

disp(datestr(night.time));

end