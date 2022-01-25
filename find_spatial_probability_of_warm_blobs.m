% This script is used to aggregate warm blob statistics sampled by RHB on different days
% and also blob statistics on different large scale condition dubbed
% (sugar-, gravel-, flowers-, fish-favored time);

% date: Jan 19, 2022
% Note: to-do: for all data;
%      - update the warm pixels (selected within the blobmaks (blobImage
%       =1), instead of using the bounding box.) 
%      - need to rethink how to remove doppelgangers. (sometimes, maybe using centroid distance instead?)
%        (or put some bounds on the maximum size of the bounding box..., sometimes, my bounding box becomes too large.)

%% load in all data:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
datasvdir = [datadir filesep 'blob_data/detrended_large_scale_SSTgrad_updated'];

blob_files = dir([datasvdir filesep '*_individual_warm_blobs_noLargeScaleSSTGrad_*.mat']);
filenames = {blob_files.name};
nfiles = length(filenames);

c0 = 1;
for i = 1:nfiles
    FN = filenames{i};
    load([datasvdir filesep FN]),'indiv_blobs_collection');
    
    nb = length(indiv_blobs_collection);
    cN = c0+nb-1;
    all_blobs(c0:cN) = indiv_blobs_collection;  
    all_blobstats(c0:cN)
    c0 = cN+1;
end

%% for all individual blobs:
% use a function to provide sorting.
grouped_blobs = group_blobs_by_U10_and_LTS_conditions(all_blobs);

% the data in above has the following fields:
% grouped_blobs.(CloudType); CloudType =
% ['Sugar','Gravel','Flowers','Fish'];

% Then I can creat a map of 2D histogram that provide the count or
% probability of having a cloud blobs in the area bin. probability map.
cloud_types = {'sugar','gravel','flowers','fish'};

% should I actually find the probability for each type?
% The probability will be higher on the clear region obviously.
for tt = 1:4
    CloudType = cloud_types{tt};
    
    blobs_tmp = grouped_blobs.(CloudType);
    blobs_stat_tmp = [blobs_tmp.stats_selected];
    blob_EqvDiam = [blobs_stat_tmp.EqvDiam];
    
    % for each cloud type: do the following:
    BlobPixelLocs = get_pixelLocs_in_boundingbox(blobs_tmp);
    %blob_centroids = vertcat(blobs_tmp.GeoLocs);
    X = BlobPixelLocs(:,1);
    Y = BlobPixelLocs(:,2);
    XBinEdge = [-59:0.25:-48];
    YBinEdge = [8:0.25:18];
    
    figure(tt);
    % this probability haven't included information on size.
    hh = histogram2(X,Y, XBinEdge, YBinEdge,'Normalization','pdf', ...
        'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
    colorbar
    title(CloudType);
    
    pdf_conditioned.(CloudType)=hh.Values';
end

% do a probability map using all data:
BlobPixelLocs_all = get_pixelLocs_in_boundingbox(all_blobs);
%blob_centroids = vertcat(blobs_tmp.GeoLocs);
X = BlobPixelLocs_all(:,1);
Y = BlobPixelLocs_all(:,2);
XBinEdge = [-59:0.25:-48];
YBinEdge = [8:0.25:18];

XBinCen = 0.5*(XBinEdge(1:end-1) + XBinEdge(2:end));
YBinCen = 0.5*(YBinEdge(1:end-1) + YBinEdge(2:end));

figure(6);
% this probability haven't included information on size.
h = histogram2(X,Y, XBinEdge, YBinEdge,'Normalization','pdf', ...
    'FaceColor','flat','DisplayStyle','tile','ShowEmptyBins','on');
colorbar
title('Probability of warm blob detection from the 2-month data');

pdf_all = h.Values;

hold on;
[C,h]=contour(XBinCen, YBinCen, pdf_blobdetection',[0.005:0.0025:0.03],'w');

save([datadir filesep 'probability_density_function_for_warmblob_detection_ATOMIC.mat'],...
    'pdf_conditioned','pdf_all','XBinCen','YBinCen','XBinEdge','YBinEdge');
