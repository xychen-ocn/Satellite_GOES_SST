`function cloudInfo = compute_cloudiness_over_features_from_L3C_cloudmasks(t, search_box, blobsIn)
% purpose: this script is used to compute cloudiness from L3C cloudmasks,
% in the region of interest. 
% Inputs: 
%     t: time of request;
%     search_box: area to compute cloudiness from L3C cloudmask; 
%  - need to deal with the fact that L3C data is about 2.5x finer than the
%  L4 data; (how do we use the cloud mask information for a coarser grid:
%  but now I don't think we really need to change the resolution.)
%
%

L3C_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/g16';
% load in L3C daily dataset;
tday = floor(t);
dateID = datestr(tday, 'mmmdd');
DOY = tday - datenum(2019, 12, 31);

% daily cloud mask.
matFN = ['GOES_SST_' num2str(DOY,'%3.3i') '_' dateID '_lon62W_TO_42W_lat5N_TO_25N.mat'];
disp(['loading' dateID])
load([L3C_datadir filesep matFN]);

[GOES_LON, GOES_LAT] = meshgrid(GOES_ATOMIC.lon(:,1), GOES_ATOMIC.lat(:,1));

SST_all = GOES_ATOMIC.sea_surface_temperature;
time_all = GOES_ATOMIC.time_num;

% check cloud mask and throw away maps that contains all cloudy pixels.
%% QC data:
% through away maps that is 100% garbage;
nt0 = size(SST_all,3); nt = nt0;
i=1;
while i<=nt
    SST_tmp = SST_all(:,:,i);
    idxs = find(isnan(SST_tmp));
    if length(idxs) == numel(SST_tmp) 
        SST_all(:,:,i) = [];
        time_all(i) = [];
        
        nt = size(SST_all,3);
    else
        i = i+1;
    end
    
end

BndBox.cen = blobsIn.GeoLocs;
BndBox.Size = blobsIn.BoundingBoxSize;  
search_RadRatio = search_box.sRadiusRatio;

for ib = 1:length(search_box.BndBox_dw)
    % get data in the search box
    if 1==0
    vertices = search_box.SBox_dw_vertices{ib};
    search_square_vx = vertices(1,:);
    search_square_vy = vertices(2,:);
    
    [in, on] = inpolygon(GOES_LON, GOES_LAT, search_square_vx, search_square_vy);
    smask = in | on;
    end
    
    %% how to do this when the search area is a circle?
     % use the center of the blob and the distance for searching.
     L = max(BndBox.Size(ib,:))*search_RadRatio;
     xb = BndBox.cen(ib,1); yb = BndBox.cen(ib,2);
     dist = sqrt((GOES_LON - xb).^2 + (GOES_LAT-yb).^2);
     smask = dist<=L*1.5;
    % extract cloud mask at selected location for the entire day:
    SST_cloudmask= [];
    for it =1:length(time_all)
        SST_tmp = SST_all(:,:,it);
        SST_cloudmask(:,it) = SST_tmp(smask);
    end
    search_lon = GOES_LON(smask); search_lat = GOES_LAT(smask);
    
    % compute cloudiness (6:3:24) for different period and see how it evolves.
    if length(time_all)>=6
        cldfreq = compute_cloudfreq(SST_cloudmask);
    else
        cldfreq = NaN;
    end
    
    % visualize:
%     if ~isnan(cldfreq)
%         BlobMask = blobsIn.blob_image{ib};
%         BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);
% 
%         
%         figure(7);
%         scatter(search_lon, search_lat, 50, cldfreq,'filled','marker','s');
%         hold on;
%         % plot the warm blob element contours on here:
%         contour(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask),'-r');
%         %hold off
%         colormap(jet);
%         colorbar;
%         pause(0.1)
%     end
    
    % save data: 
    cloudInfo(ib).SST_cloudmask = SST_cloudmask;
    cloudInfo(ib).time = time_all;
    cloudInfo(ib).cloudfreq = cldfreq;
    cloudInfo(ib).lon = search_lon;
    cloudInfo(ib).lat = search_lat;
    
    % I wanted to have the 2-month cloud frequency at the location of each
    % eddies. need to make this (2-month cloud frequency data stored for access.)
    
    
end


return