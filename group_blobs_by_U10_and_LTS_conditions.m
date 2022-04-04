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



% load ncfile;
eraFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
era.time = ncread(eraFN,'time');
era.timenum = datenum('1800-01-01') + era.time/24;
era.U10_domainMean = ncread(eraFN,'U10_domainAve');
era.LTS_domainMean = ncread(eraFN, 'LTS_domainAve');

% do daily average on two quantity
eradays = unique(floor(era.timenum));
for i = 1:length(eradays)
    tmask = era.timenum>=eradays(i) & era.timenum<eradays(i)+1;
    era.U10_dailyDomainAve(i) = mean(era.U10_domainMean(tmask));
    era.LTS_dailyDomainAve(i) = mean(era.LTS_domainMean(tmask));
end
era.day_timenum = eradays;

% use the daily average quantity to provide a condition tag;
cond_U10 = era.U10_dailyDomainAve>=8;
cond_LTS = era.LTS_dailyDomainAve>=15;

cond.flowers = cond_U10& cond_LTS;
cond.gravel = cond_U10 & ~cond_LTS;
cond.sugar = ~cond_U10 & ~cond_LTS;
cond.fish = ~cond_U10 & cond_LTS;

cloud_types = {'gravel','flowers','sugar','fish'};

for m = 1:4
    CT = cloud_types{m};
    dayType.(CT) = era.day_timenum(cond.(CT));
end
% use time as a friend:
BlobTimes = [blobsIn.time];   % hourly resolution unless missing data (etc.)

for tt = 1:4
    CloudType = cloud_types{tt};
    CN = [CloudType, '_tag'];
   
%    % this following line will not work because cond is for daily and the
%    % era timenumber is 3 hourly.   
%    idx_sel = [];
%    for id = 1:length(dayType.(CT))
%        dhere = dayType.(CT)(id);
%        tmask_1d = era.timenum>=dhere & era.timenum<dhere+1;
%        idx_sel = [idx_sel, find(tmask_1d==1)];       
%    end
%    
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