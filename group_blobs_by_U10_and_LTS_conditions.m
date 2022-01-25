function grouped_blobs = group_blobs_by_U10_and_LTS_conditions(blobsIn)
% purpose: group input blobs data by 4 types of condition using the time as
%          a way to provide grouping.
% Output: grouped blobs.

% get ERA data and type flag:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
eraFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
era.time = ncread([datadir filesep eraFN],'time');
era.timenum = datenum('1800-01-01') + era.time/24;
era.sugar_tag = ncread(eraFN,'sugar_tag');
era.fish_tag = ncread(eraFN,'fish_tag');
era.gravel_tag = ncread(eraFN,'gravel_tag');
era.flowers_tag = ncread(eraFN,'flowers_tag');


% use time as a friend:
BlobTimes = [blobsIn.time];   % hourly resolution unless missing data (etc.)

cloud_types = {'sugar','gravel','flowers','fish'};
for tt = 1:4
    CloudType = cloud_types{tt};
    CN = [CloudType, '_tag'];
    sel_time = era.timenum(era.(CN)==1);
    
    fid = fopen(['grouping_output_' upper(CloudType) '.txt'],'w');

    fprintf(fid,'%s\n', ['*****', upper(CloudType) , '*****']);
    
    % use time to sort blobs:
    idx_sel=[];
    for it = 1:length(sel_time)
        disp(['--> era at ', datestr(sel_time(it))]);
        fprintf(fid,'%s\n', ['--> era at ', datestr(sel_time(it))]);
        % if the satellite time (hrly) is within the 3hrly window of the
        % era selected time:
        %idx0 = find(time_stack==sel_time(it));
        idxs = find((BlobTimes-sel_time(it)>=0) & (BlobTimes-sel_time(it)<3/24));
        if ~isempty(idxs)
            disp(['   satellite from ', datestr(BlobTimes(idxs(1))) ...
                  ' ~ ', datestr(BlobTimes(idxs(end)))]);
              
            fprintf(fid,'%s\n', ['   satellite from ', datestr(BlobTimes(idxs(1))) ...
                  ' ~ ', datestr(BlobTimes(idxs(end)))]);
            
            idx_sel = [idx_sel, idxs];
        end
    end
    % indices used for select the right groups of blobs out
    idx_sel = unique(idx_sel);  % just in case repeated indices.  
    
    fclose(fid);
    
    grouped_blobs.(CloudType) = blobsIn(idx_sel);

end



%
return