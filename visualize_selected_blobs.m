function hfig = visualize_selected_blobs(LON, LAT, SSTin, blobsIn, cldmask, varargin)
% Purpose: this is a script that help to plot the final results of
% inidividaul SST blobs detected through varying dT thres.
% Inputs:
% Outputs:
switch nargin
    case 5
        cmapname = 'bwr';
        ncolor = 16;
    case 7
        cmapname = varargin{1};
        ncolor = varargin{2};
end
cmap = getPyPlot_cMap(cmapname,ncolor);

if max(SSTin(:))<10
    itype = 1;
else
    itype=2;
end
    

%cmap = redblue;
hfig = figure(10);
pcolor(LON, LAT, SSTin);shading flat
hb = colorbar;colormap(cmap);
%caxis([0,1]);
hold on;
if itype ==1
    [C,h] = contour(LON, LAT, SSTin,[0:0.1:0.8],'-w');
    clabel(C,h,[0:0.2:0.8],'labelspacing',600);
else
    [C,h] = contour(LON, LAT, SSTin,[25:0.5:29],'--w');
    clabel(C,h,[25:0.5:29],'labelspacing',600);

    
end

scatter(LON(cldmask), LAT(cldmask),6,[0.75 0.75 0.75],'.');
hold on;
% plot blob centers:
plot(blobsIn.GeoLocs(:,1), blobsIn.GeoLocs(:,2), '+r');
% plot bounding box for each blob center (instead of the contour for the
% connected elements)
% do no plot bounding box.
% [VX, VY] = get_boundingbox_vertices;
% plot(VX, VY, '-r')
% plot contour of the warm blobs instead:
% box_coord = get_boundingbox_coord;

for i = 1:size(blobsIn.GeoLocs,1)
    if ismember('blobImageCoord', fieldnames(blobsIn))
        contour(blobsIn.blobImageCoord(i).lon, blobsIn.blobImageCoord(i).lat, blobsIn.blobImage{i},'k');
    else
        contour(blobsIn.blob_image_GeoCoord(i).lon, blobsIn.blob_image_GeoCoord(i).lat, blobsIn.blob_image{i},'k');
    end

    % note: the box centroids seem to shift a bit. Probably because I
    % didn't ge the centroid perfectly correct?
end

xlabel('Longitude (^{\circ})');ylabel('Latitude (^{\circ}N)');
if itype==1
    caxis([-0.8, 0.8]);
    set(hb,'Ticks',[-1:0.2:1]);
else
    caxis([25.5, 28]);
    set(hb,'Ticks',[25.5:0.5:28]);
    
end
%hold off;


% nested function:
    function [VX, VY]=get_boundingbox_vertices
        
        blob_cen = blobsIn.GeoLocs;
        bbox_size = blobsIn.BoundingBoxSize;
        
        x1 = blob_cen(:,1) - 0.5*bbox_size(:,1);  y1 = blob_cen(:,2) + 0.5*bbox_size(:,2); 
        x2 = blob_cen(:,1) + 0.5*bbox_size(:,1);  y2 = y1;
        x3 = x2;  y3 = blob_cen(:,2) - 0.5*bbox_size(:,2);
        x4 = x1;  y4 = y3;
        
        VX = [x1, x2, x3, x4, x1]';
        VY = [y1, y2, y3, y4, y1]';
    end


    function coord = get_boundingbox_coord
         blob_cen = blobsIn.GeoLocs;
         bbox_size = blobsIn.BoundingBoxSize;
         
         for j = 1:length(blob_cen)
             lon0 = blob_cen(j,1)-0.5*bbox_size(j,1);
             lonN = blob_cen(j,1)+0.5*bbox_size(j,1);
             
             lat0 = (blob_cen(j,2)-0.5*bbox_size(j,2));
             latN = (blob_cen(j,2)+0.5*bbox_size(j,2));
             [ny, nx ] = size(blobsIn.blobImage{j});
             
             
             coord(j).lon = linspace(lon0, lonN, nx);
             coord(j).lat = linspace(lat0, latN, ny);
         end



%          LON_cutouts = [blobsIn.LON_cutouts];     % get all the LON location for the cutouts togeter:
%          LAT_cutouts = [blobsIn.LAT_cutouts];     % idem.
% 
%          nblobs = length(LON_cutouts);
% 
%          edge_width= 0.1;                       % this is what I did by default, cutout is 0.1° wider than the Bounding box
% 
%          %Xbox = []; Ybox=[];
%          for ib = 1:nblobs
%              % for each blob, get out the pixels location within the bounding box
%              LON_here = LON_cutouts{ib};
%              LAT_here = LAT_cutouts{ib};
%              
%              lon0 = min(LON_here(:)); lonN = max(LON_here(:));
%              lat0 = min(LAT_here(:)); latN = max(LAT_here(:));
%              
%              % build mask to get BBox pixels out:
%              lon_mask = LON_here>=(lon0+edge_width) & LON_here<=(lonN-edge_width);
%              lat_mask = LAT_here>=(lat0+edge_width) & LAT_here<=(latN-edge_width);
%              boxmask = lon_mask&lat_mask;
%              
%              % these are vectors
%              LON_BBox = LON_here(boxmask);
%              LAT_BBox = LAT_here(boxmask);
%              
%              [ny, nx] = size(blobsIn.blobImage{ib});
%              coord(ib).LON = reshape(LON_BBox, nx, ny);
%              coord(ib).LAT = reshape(LAT_BBox, nx, ny);
%              
%          end
%          
    end
end