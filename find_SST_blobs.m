%function [blob_coord,blob_image, f1, stats_sel, stats_exc] = find_SST_blobs(f, fthres, minArea, cldmask,varargin)
function [blob_out, status] = find_SST_blobs(f, fthres, minArea, cldmask,varargin)

   % ==================================================================== %
   % Purpose: this function will find SST blobs and isolate those adjacent
   % to the cloudy pixel.
   % This will get statistics about the size and shape of SST from this
   % region.
   % 
   % mostly adopted from a previous script used for identify cloud blobs.
   % Input: f [Nlat x Nlon] double
   %        f0: threshould fro f;
   %        minArea: minimum area to identify an SST blob
   %        cldmask: logical [Nlat x Nlon]
   % note: the data gaps in the input filled have been filled, but cloud
   % mask will be used to mask out any identified blob influenced by such
   % filling. (Jan 11, 2022)
   
   [Ny,Nx]=size(f);
   [XX, YY]=meshgrid([1:Nx], [1:Ny]);
  
   
   default_cloudcond_type = 'relaxed';
   default_cloudcover_thres = 0.1;
   default_checkflag=false;
   
   p = inputParser;
   % define parsing rules:
   addRequired(p,'f', @isnumeric);
   addRequired(p, 'fthres', @isnumeric);
   addRequired(p, 'minArea', @isnumeric);
   addRequired(p, 'cldmask');
   checkString = @(s) any(strcmpi(s, {'relaxed','strict'}));
   addParameter(p, 'cloudcond_type', default_cloudcond_type, checkString);
   addParameter(p, 'cloudcoverage_thres', default_cloudcover_thres, @isnumeric);
   addParameter(p, 'LON_grid', XX, @isnumeric);
   addParameter(p, 'LAT_grid', YY, @isnumeric);
   addParameter(p,'checkflag', default_checkflag, @islogical);
   
   % ready to parse the inputs:
   parse(p, f, fthres, minArea, cldmask, varargin{:});
   f = p.Results.f;
   fthres = p.Results.fthres;
   minArea = p.Results.minArea;
   cldmask = p.Results.cldmask;
   cloudcond_type = p.Results.cloudcond_type;
   cc_thres = p.Results.cloudcoverage_thres;
   XX = p.Results.LON_grid;
   YY = p.Results.LAT_grid;
   checkflag = p.Results.checkflag;
   
   
   % ---------- main code starts from here on ------------- %
  
   fo = f;      % save original data, just for ploting at the end;
   NConn =8;   %  Janssens et al. (2021) used 6-connectivity??
   
   %f=SST_anom; f0=dT_thres;
 
   BI=((f)>fthres); %matrix with ones where f excceeds threshold
   CCr = bwconncomp(BI,NConn); % defines the connected components of B
   Lr = labelmatrix(CCr); % each integer corresponds to one of the connected components
   num_events_orig = max(Lr(:));
   
   fprintf(1,[int2str(num_events_orig),' - ini number of events\n']);
   
   stats = regionprops(Lr,f,'BoundingBox','WeightedCentroid',...
       'Image','Orientation','Eccentricity','MaxFeretProperties',...
       'EquivDiameter','EulerNumber', 'Perimeter','Area', ...
       'MajorAxisLength','MinorAxisLength', 'PixelIdxList','ConvexHull','ConvexImage');           % add returning pixel indices (p-element vector)
   
   if num_events_orig >0
       
       box=[stats.BoundingBox];
       cent=[stats.WeightedCentroid]; % first elem: x-dir, second: y-dir
       orien = [stats.Orientation];
       eccen = [stats.Eccentricity];
       if ~isempty(stats)
           maxFeretAng = [stats.MaxFeretAngle];
           maxFeretCoord = reshape([stats.MaxFeretCoordinates],4,[]);
           maxFeretDiam = [stats.MaxFeretDiameter];
       end
       EqvDiam = [stats.EquivDiameter];
       EuNum = [stats.EulerNumber];
       Perim = [stats.Perimeter];
       area = [stats.Area];
       major_axislen = [stats.MajorAxisLength];
       minor_axislen = [stats.MinorAxisLength];
       
       for i = 1:num_events_orig
           pixel_idx{i} = stats(i).PixelIdxList;
       end
       
       
       box_candidates=reshape(box,4,max(Lr(:)));
       
       % ini_t=box_candidate(1,event) ; ini_x=box_candidate(2,event)
       % width_t=box_candidate(3,event);  width_x=box_candidate(4,event)
       
       cent_candidates=reshape(cent,2,max(Lr(:)));
       blob_image0 = {stats.Image};
       
       
       % in Janssens' work, the area needs to be larger than the area of 4 pixels.
       % with 4-connectivity
       cond_area = area>minArea;
       %cond_area = area>prctile(area,25);
       %disp(['-crit: area=', num2str(prctile(area,25))]);
       
       
       %%%% ----- construct a database that stores the pixel indices for edges -----
       left_edge = [1:Ny];
       top_edge = 1:Ny:(Nx-1)*Ny+1;
       right_edge = (Nx-1)*Ny:1:Nx*Ny;
       bottom_edge = Ny:Ny:Nx*Ny;
       
       edge_matrix = unique([left_edge, top_edge, right_edge, bottom_edge]);
       
       cond_edge = false(size(cond_area));
       for i = 1:num_events_orig
           % if any pixel in the originally identified object contains edge
           % pixels.
           cond_edge(i) = any(ismember(pixel_idx{i}, edge_matrix  ));       % dimension: [1x num_of_events]
       end
       
       %%%% ------ construct a database that stores the pixel indices for
       %%%% ------ cloudy-pixels in the SST map
       % 1. identify the indices for the cloudy pixels
       cloudy_indx = find(cldmask==1);
       
       % 2. for each identified object, check if the the object contains
       % cloudy pixels, if so, save the objects separately.
       cond_cloud = false(size(cond_area));
       %size(cond_cloud);
       for i = 1:num_events_orig
           if strcmpi(cloudcond_type, 'strict')
               % if an blob event contains cloudy_indx, then the blob is removed.
               cond_cloud(i) = any(ismember(pixel_idx{i}, cloudy_indx));
               
           else
               % a less strict condition:
               % display('use relaxed condition to identify cloud influenced blobs..')
               % display(['-- cloud coverage threshold = ', num2str(cc_thres)]);
               % -- check how many percentage of the blob is covered by clouds
               numpixel_blob = length(pixel_idx{i});
               numpixel_cloud = length(find(ismember(pixel_idx{i}, cloudy_indx)==1));
               
               cloud_coverage =  numpixel_cloud/numpixel_blob;
               %display(['-- cloud coverage = ', num2str(cloud_coverage)]);
               cond_cloud(i) = (cloud_coverage>=cc_thres);    % tag as cloud influenced (to be excluded.)
           end
       end
       
       %cond_lon = (box_candidates(2,:)> 1 & (box_candidates(2,:)+box_candidates(4,:))<N); % have to deal with periodicity
       
       %idx_events = find(cond_area  & ~cond_edge ); % keep only cases that satisfy the two conditions above
       idx_events = find(cond_area & ~cond_edge & ~cond_cloud);
       
       
       edge_events = find(cond_area  & cond_edge);
       cloudy_events = find(cond_area & cond_cloud);
       excluded_events= find(cond_area & (cond_cloud|cond_edge));
       
       
       if ~isempty(idx_events)
           num_events = max(size(idx_events));
       else
           num_events=0;
       end
       fprintf(1,[int2str(num_events),' -  actual number of events\n']);
       fprintf(1,[int2str(length(excluded_events)),' -  number of cloud or edge events\n']);
       
       blob_image1 = blob_image0(idx_events);
       blob_box1 = box_candidates(:,idx_events)'; %(timw lonw) . (lon, lat) 3=x, 4=y
       blob_cent1 = cent_candidates(:,idx_events)'; %(timcent loncent) (lon,lat)
       % add extra quantities here:
       stats_sel.eccen = eccen(idx_events);
       stats_sel.orientation  = orien(idx_events);
       if ~isempty(stats)
           stats_sel.maxFeretAng = maxFeretAng(idx_events);
           stats_sel.maxFeretDia = maxFeretDiam(idx_events);
           stats_sel.maxFeretCoord = maxFeretCoord(:,idx_events);
           stats_sel.MajAxisLen = major_axislen(idx_events);
           stats_sel.MinAxisLen = minor_axislen(idx_events);
       end
       stats_sel.EqvDiam = EqvDiam(idx_events);
       stats_sel.EuNum = EuNum(idx_events);
       stats_sel.Perimeter = Perim(idx_events);
       stats_sel.Area = area(idx_events);
       
       
       stats_exc.MajAxisLen= sum(area(excluded_events));
       stats_exc.MinAxisLen = EqvDiam(excluded_events);
       stats_exc.Area = area(:,excluded_events);
       stats_exc.eccen = eccen(excluded_events);
       stats_exc.orientation = orien(excluded_events);
       
       
       % I need to save the area of objects I dropped out to compute the cloud
       % density;
       
       % the selected data:
       BW_sel = ismember(Lr, idx_events);
       f1=(BW_sel>0);% for ploting only
       
       blob_coord1 = [blob_cent1 blob_box1];
       
       blob_image = blob_image1;
       blob_coord = blob_coord1;   %get better understanding of this parameter.
       
       % the excluded blobs:
       BW_exc = ismember(Lr, excluded_events);
       fexc = (BW_exc>0);
       blob_coord2 = [cent_candidates(:,excluded_events)' box_candidates(:,excluded_events)'];
       
       %blob_image2 = blob_image0(excluded_events);
       
       
  
       %hold off;
       
       
       %% convert center location from grid index to Geophysical location;
       if num_events>0
           lon_indx = 1:Nx; lat_indx=1:Ny;
           lon_vec = XX(1,:); lat_vec = YY(:,1);
           dlon = mean(diff(lon_vec)); dlat = mean(diff(lat_vec));
           
           blob_GeoLocs.clon = interp1(lon_indx, lon_vec, blob_coord(:,1));
           blob_GeoLocs.clat = interp1(lat_indx, lat_vec', blob_coord(:,2));
           excblobs.clon = interp1(lon_indx, lon_vec,blob_coord2(:,1));
           excblobs.clat = interp1(lat_indx, lat_vec,blob_coord2(:,2));
           
           % get coordinate for the blob image: (and convert it to the
           % geolocations);
           bbox_pixel_coord = build_coordinate_for_blobImage(blob_box1);
           bbox_GeoCoord = pixelCoord_to_geoCoord(bbox_pixel_coord);
           
           
           % PLOT CHECK
           if checkflag
               figure(10); clf;
               pcolor(XX, YY, fo);shading flat%caxis([0,1]);
               colorbar;colormap('parula');
               hold on;
               contour(XX,YY, fo,[-0.02,0.02],'-w');
               scatter(XX(cldmask), YY(cldmask),6,[0.75 0.75 0.75],'+');
               hold on;
               if num_events>0
                   contour(XX, YY, f1,'r');
                   % construct the inidividual blob contour (to doboule check)
                   % tested, it can give the correct warmer blob pixels.
                   %             for ib = 1:length(blob_image)
                   %                 pcolor(bbox_GeoCoord(ib).lon, bbox_GeoCoord(ib).lat, double(blob_image{ib}));shading flat;
                   %             end
                   
                   
                   hold on;
                   plot(blob_GeoLocs.clon,blob_GeoLocs.clat,'r+','markersize',10);
                   
                   contour(XX, YY,fexc,'--k','linewidth',1.);
                   plot(excblobs.clon, excblobs.clat, '+k','markersize',8);
               end
               
               
               xlabel('lon (# grid points)');ylabel('lat(# grid points)');
               title('f, blobs, and centroids')
               hold off;
           end
           %% store all output in 1 single structure.
           blob_out.GeoLocs = [blob_GeoLocs.clon, blob_GeoLocs.clat];          % first column: lon, second column: lat
           blob_out.BoundingBoxSize = [blob_coord(:,5)*dlon , blob_coord(:,6)*dlat];
           blob_out.stats_selected = stats_sel;
           %blob_out.stats_excluded = stats_exc;  # I don't need this  infomation yet.
           blob_out.blob_mask = f1;
           blob_out.blob_image = blob_image;                                   % N cells
           blob_out.blob_image_GeoCoord = bbox_GeoCoord;
           blob_out.PixelIds = pixel_idx;

           status = 1;
       else
           % no cloud events
           status = -1;
       blob_out= NaN;
       end
           
   else
       % no init events.
       status = -1;
       blob_out= NaN;
       blob_coord = [];
       blob_image=[];
       f1 = []; 
       stats_sel=[]; 
       stats_exc=[];
        return
       
       
   end
        
   
   %% nested function
    function pixel_coord = build_coordinate_for_blobImage(blob_box)
        
        x_lowleft = blob_box(:,1);
        y_topleft = blob_box(:,2);
        width = blob_box(:,3);
        height = blob_box(:,4);
        
        nb = length(x_lowleft);
        for k = 1:nb
            pixel_coord(k).x = round(x_lowleft(k)):round(x_lowleft(k))+width(k)-1;
            pixel_coord(k).y = round(y_topleft(k)):round(y_topleft(k))+height(k)-1;
            
        end
       
    end

    function GeoCoord = pixelCoord_to_geoCoord(blob_pixel_coord)
        %lon_indx = 1:Nx; lat_indx=Ny:-1:1;
        %lon_vec = XX(1,:); lat_vec = YY(:,1);
        %dlon = mean(diff(lon_vec)); dlat = mean(diff(lat_vec));
        
        nb = length(blob_pixel_coord);
        for k =1:nb
            GeoCoord(k).lon = interp1(lon_indx, lon_vec, blob_pixel_coord(k).x);
            GeoCoord(k).lat = interp1(lat_indx, lat_vec, blob_pixel_coord(k).y);
        end
        
    end
        
end