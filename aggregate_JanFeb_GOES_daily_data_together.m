% save all the daily g16 data together in one dataset:
clear all; clc; 

datadir = 'g16';
files = dir([datadir filesep 'GOES_SST_*.mat']);
filenames = {files.name};

it0=1;
for i = 1:length(files)
    load([datadir filesep filenames{i}]);
    
    datafields = fieldnames(GOES_ATOMIC);
    
    GOES_data.lon = GOES_ATOMIC.lon(:,1);
    GOES_data.lat = GOES_ATOMIC.lat(:,1);
    
    nt = length(GOES_ATOMIC.time);
    itN = it0+nt-1;
    disp(['it0=' num2str(it0), '; itN=' num2str(itN) '; nt=', num2str(nt)]);
    if nt~=24
    pause
    end
    
    % stack data together for all variables:
    for j = 1:length(datafields)
        FN = datafields{j};
        if length(size(GOES_ATOMIC.(FN)))==3 

            GOES_data.(FN)(:,:,it0:itN) = GOES_ATOMIC.(FN);
            
        elseif length(size(GOES_ATOMIC.(FN)))==2
            if numel(GOES_ATOMIC.(FN)) == length(GOES_ATOMIC.(FN))
            GOES_data.(FN)(it0:itN) = GOES_ATOMIC.(FN);
            end
        end
    end
    
    it0 = itN+1;
end

% save data:
save([datadir filesep 'GOES_SST_JanFeb2022_ATOMIC_region.mat'],'GOES_data','-v7.3');
