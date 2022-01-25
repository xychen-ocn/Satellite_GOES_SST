function [VX, VY]=get_boundingbox_vertices(blobsIn)
% Input: blobsIn (individual blobs outputted from remove_doppelgangers)

blob_cen = blobsIn.GeoLocs;
bbox_size = blobsIn.BoundingBoxSize;

x1 = blob_cen(:,1) - 0.5*bbox_size(:,1);  y1 = blob_cen(:,2) + 0.5*bbox_size(:,2);
x2 = blob_cen(:,1) + 0.5*bbox_size(:,1);  y2 = y1;
x3 = x2;  y3 = blob_cen(:,2) - 0.5*bbox_size(:,2);
x4 = x1;  y4 = y3;

VX = [x1, x2, x3, x4, x1]';
VY = [y1, y2, y3, y4, y1]';
end