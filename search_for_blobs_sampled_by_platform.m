function [PFM_sampled_blobs, obsv_data, status] = search_for_blobs_sampled_by_platform(platform_name, blobs)
% Purpose: 1. fetch data by platform name
%          2. using the platform data segments to locate the warm blobs
%          sampled by the platform during ATOMIC.
% Inputs: platform_name (string)
%         blobs (structure, output from remove_doppelgangers)
% Outputs: PFM_sampled_blobs (warm blobs that the platform trajectorys
%          crossed.
% Date   : Jan 17, 2022 (drafted)

%% load in data stored in the following directory;
ATOMIC_dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/ATOMIC/rhb/object_oriented_scripts';
addpath(ATOMIC_dataroot);
datadir = [ATOMIC_dataroot filesep 'processed_data'];
dataFN = [platform_name '_along_wind_segments_ready_for_wavelet_coherence.mat'];


load([datadir filesep dataFN]);
eval(['segobj =', platform_name '_alongwind_segobj;']);     % change variable name to a universal name.


% take out platform data according to blob time:
obsv_data = take_out_platform_data_according_to_blobTime(segobj, blobs);

%% search for blobs in which the platform trajectory passed through;
%  search for blobs in specific time window (1 hr) becuase the blob is
%  updated hourly. (t_blob-0.5: t_blob+0.5) (center on t_blob)
%  It seems easier to refine search region if I consider specific time
%  window. 

% 1. narrow down the seach region:
% find blobs with centroid located in a 1 degree radius from the
% observation track
blob_cen = blobs.GeoLocs;
tmask = ((obsv_data.time-blobs.time)*24<=0.5) & ((obsv_data.time-blobs.time)*24>=-0.5);            % 
obsv_avelon = mean(obsv_data.lon(tmask),'omitnan');
obsv_avelat = mean(obsv_data.lat(tmask),'omitnan');
search_radius = 1;    % degree

x1 = obsv_avelon-search_radius; y1 = obsv_avelat+search_radius;
x2 = obsv_avelon+search_radius; y2 = y1;
x3 = x2;                        y3 = obsv_avelat-search_radius;
x4 = x1;                        y4 = y3;

search_region_XV = [x1, x2, x3, x4];
search_region_YV = [y1, y2, y3, y4];
%figure
%patch(search_region_XV, search_region_YV, 'r','FaceColor','none');
%hold on;
%plot(obsv_data.lon(tmask), obsv_data.lat(tmask),'.-k');

blobInFlag = inpolygon(blob_cen(:,1), blob_cen(:,2), search_region_XV, search_region_YV);

% I need a utility function to subset blob information 
blobInIndx = find(blobInFlag==1);
nblobs_in = length(blobInIndx);
disp([num2str(nblobs_in) ' blobs are found in the search radius.']);

if nblobs_in>=1
    selected_blobs = subset_blobs_with_indices(blobs, blobInIndx); %%
else
    status = false;
    PFM_sampled_blobs = NaN;
    return
end

%plot(selected_blobs.GeoLocs(:,1), selected_blobs.GeoLocs(:,2),'+r');

% 2. in those blob candidates, search for the ones that the platform passes by.
[XV, YV] = get_boundingbox_vertices(selected_blobs);
%plot(XV,YV,'--r');



bid = [];
for i = 1:nblobs_in
    % need to get the indx that can get the correct warm blobs.
    % a section of the RHB track (from 1 hour before the blob sampled time
    % to the blob sampled time);
    inflag = inpolygon(obsv_data.lon(tmask), obsv_data.lat(tmask), XV(1:4,i), YV(1:4,i));
    
    if any(inflag)   % there is at least 1 obsv trajectory point is in the blob bounding box.
        % count this warm blob
        bid = [bid, i];    % use this index to get blobs out 
    end
    
end

status = ~isempty(bid);       % found blobs.
if status
    disp([num2str(length(bid)) ' blob(s) found!']);
    if nblobs_in > length(bid)
        PFM_sampled_blobs = subset_blobs_with_indices(selected_blobs, bid);
        %PFM_sampled_blobs.time = blobs.time;
        
    else
        PFM_sampled_blobs = selected_blobs;
    end
else
    PFM_sampled_blobs = NaN;
    disp('no blobs near observation.');
end
%PFM_sampled_blobs.GeoLocs = selected_blobs.GeoLocs(bid,:);



end


%%%%%%%%%%%%%%%%%%%%%%% nested function 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obsv_data = take_out_platform_data_according_to_blobTime(segobj, blobs)
nsegs = length(segobj);
time_start = nan(nsegs, 2);
time_end = nan(nsegs, 2);
for iseg = 1:nsegs
    nsections = length(segobj(iseg).Value);
    for k = 1:nsections
        time_start(iseg, k) = segobj(iseg).Value(k).time(1);
        time_end(iseg, k) = segobj(iseg).Value(k).time(end);
    end
end

% use the information to get the right segment according to the day of
% interest.
blob_date = datestr(blobs.time,'mmmdd');
disp(['blob sampled time is' blob_date]);
% check if the blob time is within the start and end time of a segment?
cond = (blobs.time - time_start) .* (blobs.time - time_end);
within_flag = cond<=0;
idx = find(within_flag==1);
% convert index back to subscript to retrive the right segment;
[seg_key, sect_key]=ind2sub(size(within_flag),idx);
obsv_data = segobj(seg_key).Value(sect_key);
end

%%%%%%%%%%%%%%%%%%%%%%% nested function 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function subset = subset_blobs_with_indices(blob_struc, ind)

fieldN_list = fieldnames(blob_struc);
for i = 1:length(fieldN_list)
    FN = fieldN_list{i};
    
    if isnumeric(blob_struc.(FN))
        
        if numel(blob_struc.(FN)) > length(blob_struc.(FN))   % arrays:
            subset.(FN) = blob_struc.(FN)(ind,:);
        else % vectors or a single number
            if length(blob_struc.(FN))~=1
                if iscolumn(blob_struc.(FN)) 
                    subset.(FN) = blob_struc.(FN)(ind,:);
                else % row vec
                    subset.(FN) = blob_struc.(FN)(:,ind);
                end
            else
                subset.(FN) = blob_struc.(FN);
            end
        end
        
    elseif iscell(blob_struc.(FN))         % know apriori that my cell arrays are 1xnblobs.
        for k =1:length(ind)
            
            subset.(FN){k} = blob_struc.(FN){ind(k)};
        end
        
    elseif isstruct( blob_struc.(FN) )
        if strcmp(FN, 'stats_selected')
            stats_vars = fieldnames(blob_struc.(FN));
            for s = 1:length(stats_vars)
                SN = stats_vars{s};
                % gather up statistics associated with all the blobs:
                % after "remove_doppelgangers.m", these statistic variables are
                % 1 dimensional.
                subset.(FN).(SN) = blob_struc.(FN).(SN)(ind);
                
            end
        elseif strcmp(FN, 'blobImageCoord')
           subset.(FN) = blob_struc.(FN)(ind);
        end
           
    end     % end variable types.
end

% add other variables (without typing them out.)
% check field variables:
% subset.GeoLocs = blob_struc.GeoLocs(ind,:);
% subset.BoundingBoxSize = blob_struc.BoundingBoxSize(ind,:);



end