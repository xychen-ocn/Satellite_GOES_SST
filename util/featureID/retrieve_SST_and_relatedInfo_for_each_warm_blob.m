function blobs =retrieve_SST_and_relatedInfo_for_each_warm_blob(blobs, LON, LAT, SST, SST_LS_trend, varargin)
% purpose: cut out the SST map according to Bounding Boxes for warm blobs
%          compute SST anomaly for the cutout region as well
%          associate the blob with its large scale background SST and
%          maximum SST anomaly.
% 
% Inputs: blobs (structure output with field constructed from
%                 "remove_doppelgangers.m", this contains all the
%                 individual blobs of interest at given time instance)
%         LON, LAT (2D grid location)
%         SST  (SST map)
%         blobmask (from bounding box.)
%         
% Outputs:
% Authors:
% Date: drafted Jan 11, 2022.
%       updated Apr 6, 2022: modify cutout_extend to be certain times of
%       the equiv. radius. 
%
%SSTtype = 'warm';
switch nargin
    case(5)
        cutout_extend = 1;       % 1x EqvRad.
    case(6)
        cutout_extend = varargin{1};
%     case(7)
%         cutout_extend = varargin{1};
%         SSTtype = varargin{2};
end

NBlobs = length(blobs.GeoLocs);


% search for individual warm blob:
for i = 1:NBlobs
    
    clon = blobs.GeoLocs(i,1);   % top-left corner of the box
    clat =  blobs.GeoLocs(i,2);   % top-left corner of the box
    LX = blobs.BoundingBoxSize(i,1);
    LY = blobs.BoundingBoxSize(i,2);

    EqvDiam_dgr = blobs.stats_selected.EqvDiam(i) * blobs.L4_xres/111E3;   % approximately correct. (It is okay because I do not need super accurate here.)
    ExtRad = EqvDiam_dgr/2*cutout_extend ;

    lon0 = clon - LX/2; %- cutout_extend;   % index of the end of the box.
    lonN = clon + LX/2; %+ cutout_extend;
    
    lat0 = clat - LY/2; %- cutout_extend;
    latN = clat + LY/2; %+ cutout_extend;
    
    % create mask:
    cutout_mask = false(size(SST));
    cutout_mask((LON>(lon0-ExtRad))&(LON<(lonN+ExtRad)))= true;
    cutout_mask((LAT<(lat0-ExtRad))|(LAT>(latN+ExtRad))) = false;
    
    bbox_mask = false(size(SST));
    bbox_mask((LON>lon0)&(LON<lonN)) = true;
    bbox_mask((LAT<lat0)|(LAT>latN)) = false;
    
    % use the index to slice the SST matrix:
    SSTcutout = SST(cutout_mask);
    if ~isempty(SSTcutout)
    %SSTacutout = SST_anom(cutout_mask);
    LONcutout = LON(cutout_mask);
    LATcutout = LAT(cutout_mask);
    SST_LSG_cutout = SST_LS_trend(cutout_mask);                            % LSG = Large Scale Gradient.
    %blobmask_cutout = blobmask(cutout_mask);                                % this is not needed now, since it doesn't really exist.
    SST_LSG_bbox = SST_LS_trend(bbox_mask);
    
    latslice = LAT(:,1);
    nrow = length(find((latslice>=lat0-ExtRad)&(latslice<=latN+ExtRad)));
   
    
    blobs.LON_cutouts{i} = reshape(LONcutout, nrow, []);
    blobs.LAT_cutouts{i} = reshape(LATcutout, nrow, []);
    blobs.SST_cutouts{i} = reshape(SSTcutout, nrow, []);
    %blobs.SSTa_cutouts{i} = reshape(SSTacutout, nrow, []);
    
    
    % compute the SST anomaly in the cutout regions:
    SST_background = reshape(SST_LSG_cutout, nrow, []);
    SST_anom = blobs.SST_cutouts{i} - SST_background;
    blobs.SSTa_cutouts{i} = SST_anom;
    
    % a more general one:
%     max_SSTanom = max(SST_anom(:));
%     min_SSTanom = min(SST_anom(:));
    % 03/22/2022 updated: 4/5 updated again.
    % use the SST anom at the object centroid to decide. 
    if ~any(isnan(blobs.GeoLocs(i,:)))   % seems very strange. why??
        SST_anom_at_bcen = interp2(blobs.LON_cutouts{i},  blobs.LAT_cutouts{i}, SST_anom, ...
            blobs.GeoLocs(i,1), blobs.GeoLocs(i,2));
        if SST_anom_at_bcen> 0
            % warm anomaly:
            blobs.max_SSTa(i) = max(SST_anom(:));    % max_SSTa: maixmum magnitude of SSTa.
        elseif SST_anom_at_bcen < 0
            % cold anomaly:
            blobs.max_SSTa(i) = min(SST_anom(:));
        end
        
        % put in background SST (averaged in the bounding box region)
        nrow2 = length(find((latslice>=lat0)&(latslice<=latN)));
        SST_background_bbox = reshape(SST_LSG_bbox, nrow2, []);
        SST_bg_ave = mean(SST_background_bbox(:),'omitnan');
        blobs.ave_SSTbg(i) = SST_bg_ave;
        
    else
        blobs.max_SSTa(i) = NaN;
        blobs.ave_SSTbg(i) = NaN;
    end
    
    
    else
        blobs.LON_cutouts{i} =NaN;
        blobs.LAT_cutouts{i} = NaN;
        blobs.SST_cutouts{i} = NaN;
        blobs.SSTa_cutouts{i} = NaN;
        
        blobs.max_SSTa(i) = NaN;
        blobs.ave_SSTbg(i) = NaN;
        
    end
    
    % double check:
%     if strcmp(SSTtype,'warm')
%         if max(SST_anom(:))<0
%             disp('error: SST_anom should not be negative! check code!')
%             return
%         else
%             blobs.max_SSTa(i) = max(SST_anom(:));                          % This should always be positive isn't it??
%         end
%         
%     else
%         blobs.max_coldSSTa(i) = min(SST_anom(:));
%         
%     end
    

    
   
    
    % plot it out to confirm: (confirmed during initial development)
%     if 1==0
%         figure;
%         pcolor(blobs.LON{i}, blobs.LAT{i}, blobs.SST{i});
%         shading flat;
%         colorbar;
%         hold on
%         plot(clon, clat,'+r');
%         contour(blobs.LON{i}, blobs.LAT{i}, ...
%             reshape(blobmask_cutout,nrow,[]),'-k');
%         caxis([300-1.5, 300+1]);
%         axis('equal');
%         pause
%     end
    
end

return