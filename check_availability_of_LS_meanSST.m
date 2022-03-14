function [avail_flag,meanSST, meanSST_LS] = check_availability_of_LS_meanSST(SST_mean, times_q)
% Purpose: this is a untility function intend to check if LS temporal mean
%          SST exist for use during the input time (hourly);
% Inputs:
%        SST_LS_times, times_q, aveMethod
%
% Outputs:
%        avail_flag, LSid
%
% Date: Jan 30, 2022 (XYC)
% Note: 
%   for 3, 5, 7-day moving average, compare if the day of the request times are
%   found in the SST_LS_times;
%
%   for 15, 30-day (fixed window) average, check if the request times is in
%   any of the time groups that the averaging is performed with.
%
aveMethod = SST_mean.aveMethod;
req_tday = unique(floor(times_q));     % simply one day. implicit assumption, automatically satisfied by my code.

if strcmpi(aveMethod, 'moving')
    avail_days = floor(SST_mean.time);
    
    avail_flag=ismember(req_tday, avail_days);
    
    if avail_flag
        % find the specific index for the requested day:
        LSid = find(avail_days==req_tday);
        meanSST_LS = SST_mean.LPF(:,:,LSid);
        meanSST = SST_mean.values(:,:,LSid);
    else
        meanSST_LS = NaN;
        meanSST = NaN;
    end
    
elseif strcmpi(aveMethod, 'fixed')
    avail_days = floor(SST_mean.time);
    avail_days(isnan(avail_days))=0;
    
    avail_flag = ismember(req_tday, avail_days);
    
    if avail_flag
        % find the group the req_tday is in:
        LSid = find(avail_days==req_tday);
        [irow, icol]=ind2sub(size(SST_mean.time), LSid);
        
        meanSST_LS = SST_mean.LP(:,:,icol);
        meanSST = SST_mean.values(:,:,icol);
    else
        meanSST_LS = NaN;
        meanSST = NaN;
    end
    
end
    
    
return