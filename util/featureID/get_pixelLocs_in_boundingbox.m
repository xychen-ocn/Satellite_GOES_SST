function PixelLocs = get_pixelLocs_in_boundingbox(blobsIn)
% purpose: get the pixel location within the warm blobs (approximated by
% the bounding box for now) for probability calculation.

 LON_cutouts = [blobsIn.LON_cutouts];     % get all the LON location for the cutouts togeter:
 LAT_cutouts = [blobsIn.LAT_cutouts];     % idem.
 
 nblobs = length(LON_cutouts);
 
 edge_width= 0.1;                       % this is what I did by default, cutout is 0.1° wider than the Bounding box
 
 Xbox = []; Ybox=[];
 for ib = 1:nblobs
     % for each blob, get out the pixels location within the bounding box
     LON_here = LON_cutouts{ib};
     LAT_here = LAT_cutouts{ib};
     
     lon0 = min(LON_here(:)); lonN = max(LON_here(:));
     lat0 = min(LAT_here(:)); latN = max(LAT_here(:));
     
     % build mask to get BBox pixels out:
     lon_mask = LON_here>=(lon0+edge_width) & LON_here<=(lonN-edge_width);
     lat_mask = LAT_here>=(lat0+edge_width) & LAT_here<=(latN-edge_width);
     boxmask = lon_mask&lat_mask;
     
     % these are vectors
     LON_BBox = LON_here(boxmask);
     LAT_BBox = LAT_here(boxmask);
     
     % sanity check;
%      figure(5);
%      plot(LON_here, LAT_here,'.k');
%      hold on;
%      plot(LON_BBox, LAT_BBox, '.r');
     %hold off
     
    % pause;
     Xbox = [Xbox; LON_BBox];
     Ybox = [Ybox; LAT_BBox];
     
 end
 
 PixelLocs = [Xbox, Ybox];


return
