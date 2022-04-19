function selids = select_blobs_within_roi(blobsIn, search_region)
% Purpose: select features within region of interest;
% Inputs: 
% Outputs: subset of input blobs;

x1 = search_region(1,1);
x2 = search_region(1,2);
y1 = search_region(2,1);
y2 = search_region(2,2);

svx = [x1, x1, x2, x2, x1];
svy = [y1, y2, y2, y1, y1];

blobProps = [blobsIn.properties];
allblob_locs= vertcat(blobProps.GeoLocs);

inmask = inpolygon(allblob_locs(:,1), allblob_locs(:,2), svx, svy);
selids = find(inmask==1);    % use these indices to get the subsets of blobs;

% now, how to store the data???
% blobsIn_sub.properties = blobProps(selids);
% 




return