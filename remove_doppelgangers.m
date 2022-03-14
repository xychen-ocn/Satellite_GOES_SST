function [blobs_sel, Nblobs]=remove_doppelgangers(blobs, varargin)
% Description: This function is used to remove the multiple captures of the
% same blobs due to a change in the dT threshold. The removal is based on
% distance between the centroid of each input warm blobs and its bounding
% box size.
% Input: blobs (structure of length nthres), the field in this structure is
%        output by find_SST_blobs.m (the input is obtained at the same time
%        instance)
% Output: blobs_sel (the individual blobs selected), 
%         Nblobs: number of individual blobs.
% Author: X.Chen
% Date  : drafted on Jan 17, 2022
%         Jan 26, 2022 (update to make sure that the blob elements doesn't
%         over lap.
global LON LAT

switch nargin
    case 1
        debug=false;             % debug is off by default.
    case 2
        debug = varargin{1};      % input to turn it on.
end
blob_cenlon_tmp = []; blob_cenlat_tmp=[];
LXs_tmp = []; LYs_tmp = [];

nthres = length(blobs);
maxBoxEqvDiam_limiter = 2.0;                       % units: dgr (~200km)

%% 1. aggregate all the blobs caputred by n dT thresolds:
    % a. for given dT_thres i, compare the center of the warm blobs detected.
    % see if they form a cluster cloud. if so, the center of the cluster is
    % defined as the mean location of the warm blobs, at the given threshould.
    blob_center = vertcat(blobs.GeoLocs);                              % this is now a nblob_tot x 2 matrix
    blob_BBoxSize  = vertcat(blobs.BoundingBoxSize);                       % this matrix is also nblob_tot x 2matrix [LXs, LYs]
    blob_PixelIds_all = [blobs.PixelIds];
    blob_Images =  [blobs.blob_image];
    blobImageCoord = [blobs.blob_image_GeoCoord];
    
    blob_cenlon_tmp = blob_center(:,1);
    blob_cenlat_tmp = blob_center(:,2);
    
    LXs_tmp = blob_BBoxSize(:,1);              % units: dgr
    LYs_tmp = blob_BBoxSize(:,2);
    
    % b. for given time, the warm blob center can changes too but hopefully not
    % too much, this will allow us to determine the size and shape of the warm
    % blob by using the largest size and shape.
    cluster_tag = clusterdata(blob_center, 'Linkage','centroid','MaxClust',10);

    if debug
        figure(1);clf;
        % draw the SST cutout when dT_thres = 0;
        % sort by size:
        subplot(2,1,1);
        
        plot(blob_cenlon_tmp, blob_cenlat_tmp,'.k');
        hold on;
        % plot the blob container to outline the region occupied by the blobs.
        x1 = blob_cenlon_tmp - 0.5*LXs_tmp;  y1 = blob_cenlat_tmp + 0.5*LYs_tmp;
        x2 = blob_cenlon_tmp + 0.5*LXs_tmp;  y2 = y1;
        x3 = x2;  y3 = blob_cenlat_tmp - 0.5*LYs_tmp;
        x4 = x1;  y4 = y3;
        
        bbox_vertices_x = [x1, x2, x3, x4]';
        bbox_vertices_y = [y1, y2, y3, y4]';
        patch(bbox_vertices_x, bbox_vertices_y, 'r','FaceColor','none');
        %axis('equal')
        title('detected blobs with varies \DeltaT thresholds');
        colorbar
        
        % define regions to remove the same blobs
        %cenloc=[blob_cenlon_tmp, blob_cenlat_tmp];  % This is basically blob center
        
        subplot(2,1,2);
        scatter(blob_center(:,1), blob_center(:,2), 20, cluster_tag,'filled');
        colorbar
        colormap(jet);
        hold on;
        patch(bbox_vertices_x, bbox_vertices_y, 'r','FaceColor','none');
       % axis('equal')
        title('clustered region to remove warm blobs doppelgangers');
    end
        
        %% ready to remove DPGangers in each region:
        Ldiag = sqrt((0.5*LXs_tmp).^2 + (0.5*LYs_tmp).^2);
        nclust =max(cluster_tag);
        blobIndices=[];
        for ireg = 1:nclust
            % in each cluster: check distance between the centroid of the warm
            % blobs;
            % criterion of DPGers: distance between centroid is less than
            % 90% of sqrt((0.5LX)^2 + (0.5LY)^2);
            indv_blobIndx=[];
            % 1. get data in the cluster:
            Pntmask = cluster_tag == ireg;
            PntIdx = find(Pntmask);
            cluster_pnts0 = blob_center(Pntmask,:);
            BboxSize = Ldiag(Pntmask);
            BoxHfLx = 0.5*LXs_tmp(Pntmask);                                % half width of the bounding boxes in the current cluster 
            BoxHfLy = 0.5*LYs_tmp(Pntmask);
            
            for i = 1:length(PntIdx)
                blobPixels{i} = blob_PixelIds_all{PntIdx(i)};
                blobMasks{i} = blob_Images{PntIdx(i)};
            end
            blobMaskCoord = blobImageCoord(PntIdx);

            clear i
            % 2. compute distance between warm blob centroids: check from the
            % centroid with the largest bounding box: compute distance of the rest
            % of the warm blob centroids in the cluster to it.
            
            %% Jan 26: put a limiter on the size of the warm blobs; 2° diameter (200km)
            %  if the equivalent diameter of the bounding box
            %  (sqrt(4xArea_bbox/Pi)) > 200km --> then, remove the largest
            %  box instead and keep all the smaller elements within it.
            %
            
            disp(['input num. of blobs=' num2str(length(BoxHfLx))]);
            search_flag=true;
            cluster_pnts = cluster_pnts0;
            num_indv_wblob0 = length(cluster_pnts);
            isearch = 0;
            indices =1:length(BboxSize);
            
            while search_flag
                % search through all the blobs in each cluster.
                isearch=isearch+1;
                [maxL, maxId] = max(BboxSize);                             % BboxSize matrix will be masked, only keep the blobs indicated by ind.
                deltaX = abs(cluster_pnts(:,1) -cluster_pnts(maxId,1));
                deltaY = abs(cluster_pnts(:,2) -cluster_pnts(maxId,2));
                dist = sqrt(deltaX.^2 + deltaY.^2);
                
                maxBoxHfLx = BoxHfLx(maxId);
                maxBoxHfLy = BoxHfLy(maxId);
                
                %ThisBlobPixels = blobPixels{maxId};
                ThisBlobMask =blobMasks{maxId};
                ThisBlobCoord = blobMaskCoord(maxId);
                lonmask = LON>=ThisBlobCoord.lon(1)&LON<=ThisBlobCoord.lon(end);
                latmask = LAT>=ThisBlobCoord.lat(1)&LAT<=ThisBlobCoord.lat(end);
                bmask = lonmask&latmask;
      
                % retrieve indice: 
                ThisBlobPixelIds = find(bmask==1);
                ThisBlobElmIds = ThisBlobPixelIds(ThisBlobMask==1);

                if debug
                    figure(1); hold on; subplot(2,1,2); hold on;
                    plot(LON(bmask), LAT(bmask),'.g');
                end
                
                maxEqvDiam = sqrt(4*(2*maxBoxHfLx * 2*maxBoxHfLy)/pi);
                
                if maxEqvDiam > maxBoxEqvDiam_limiter
                    disp(['maxEqvDiam = ' num2str(maxEqvDiam, '%4.1f') ' >limiter (2dgr)']);
                    
                    % remove this big blob and keep the smaller ones;
                    rm_ind = maxId;                    % indices to be removed.
                    if isearch ==1
                        ind = indices;                     %
                    %else 
                        % the ind variable exists already from the previous
                        % search. 
                        % ind =;
                    end
                    % remove the maxId stored in the ind matrix;
                    ind(ind==maxId) = [];                   % indices to keep during this search (everything else other than the biggest box);
                    
                    num_indv_wblob = (isearch-1) + length(ind);            % do not count the current largest box.
                    
                    %disp(['num individual blob=', num2str(num_indv_wblob)]);
                    
                else
                    %disp(['maxEqvDiam = ' num2str(maxEqvDiam, '%4.2f')]);
                    cond_Lx = deltaX>maxBoxHfLx; cond_Ly = deltaY>maxBoxHfLy;
                    cond_Ldiag = dist>maxL;
                    
                    % the blob centroids are outside of the biggest (selected) bounding box.
                    ind = find(cond_Lx|cond_Ly|cond_Ldiag);     % these are the blobs that are considered different from the maximum blob.
                    rm_ind = find(~(cond_Lx|cond_Ly|cond_Ldiag));    % these are the doppelgangers, included the current maxPntIdx.
                    
                    %% add stricter check for bounding box:
                    
                    % check if there are overlapping between bounding boxes (selected by ind):
                    %  -> if so, further check if the actual blob elements overlapped.
                    %    --> if so, remove the smaller blob and keep the large ones.
                    %
                    % a. compute the distance between the current biggest bounding
                    %    box and the rest of the ones whose centroids are outside of the current bbox.
                    
                    %whos distpair_thres deltaY
                    %for ii = 1:length(ind)
                    loop = true;
                    ii=0;
                    while loop && length(ind)~=0
                        ii=ii+1;
                        % how did this error happen? Note: ind changed size
                        % when overlapped elements are found.
                        distpair_xthres=maxBoxHfLx + BoxHfLx(ind(ii));
                        distpair_ythres=maxBoxHfLy + BoxHfLy(ind(ii));
                        
                        if ~(deltaX(ind(ii))>=distpair_xthres || deltaY(ind(ii))>=distpair_ythres)
                            % overlapping b.box.
                            disp('! bounding box overlapped')
                            % check if blob elements overalpped: using the
                            % ismember function.
                            BlobAMask =blobMasks{ind(ii)};
                            BlobACoord = blobMaskCoord(ind(ii));
                            lonmask = LON>=BlobACoord.lon(1)&LON<=BlobACoord.lon(end);
                            latmask = LAT>=BlobACoord.lat(1)&LAT<=BlobACoord.lat(end);
                            bmaskA = lonmask&latmask;
                            % retrieve indice:
                            BlobAPixelIds = find(bmaskA==1);
                            BlobAElmIds = BlobAPixelIds(BlobAMask==1);
                            % this is not quite right: (the pixelIds is not
                            % right somehow.
                            overlapped = any(ismember(BlobAElmIds, ThisBlobElmIds));     % any pixels in blob A is also in this large blob under check:
                            
                            if debug
                                figure(1);hold on
                                subplot(2,1,2); hold on;
                                plot(LON(BlobAElmIds), LAT(BlobAElmIds),'.k');
                                patch(bbox_vertices_x(:,PntIdx(ind(ii))), bbox_vertices_y(:,PntIdx(ind(ii))),'r', 'FaceColor','None');
                                hold on
                                plot(LON(ThisBlobElmIds), LAT(ThisBlobElmIds),'.r');
                                pause(0.01)
                            end
                            
                            if overlapped
                                disp('! blob elements overlapped');
                                % remove the smaller blob:
                                rm_ind(end+1) = ind(ii);
                                ind(ii)=[];
                                %else
                                % nothing needs to be done.
                            end
                            
                        end
                        
                        % the ind length can change;
                        if ii>=length(ind)
                           loop = false;
                        end
                    end
                    clear ii
                    
                    
                    num_indv_wblob = (isearch-1) + 1 + length(ind);        % 1: is counting the largest warm blob, (isearch - 1) is counting the largest warm blobs from the previous searches.
                    %disp(['num individual blob=', num2str(num_indv_wblob)]);
                    indv_blobIndx=[indv_blobIndx, PntIdx(maxId)];          % keep the maximum box index.
                    
                end        % conditional fork based on Bounding box size.
                
                
                %% summarize search results:
                if debug
                    figure(1);hold on
                    subplot(2,1,2); hold on;
                    patch(bbox_vertices_x(:,PntIdx(maxId)), bbox_vertices_y(:,PntIdx(maxId)),'r', 'FaceColor','None');
                    patch(bbox_vertices_x(:,PntIdx(ind)), bbox_vertices_y(:,PntIdx(ind)),'b', 'FaceColor','None');
                    %pause
                    for ii = 1:length(ind)
                        %plot(LON(blobPixels{ind(ii)}), LAT(blobPixels{ind(ii)}),'.y');
                        % check if the Image coordinate works better:
                        % if so, use that to back out the pixel indices for
                        % overlapping check.
                        BlobACoord = blobMaskCoord(ind(ii));
                        lonmask = LON>=BlobACoord.lon(1)&LON<=BlobACoord.lon(end);
                        latmask = LAT>=BlobACoord.lat(1)&LAT<=BlobACoord.lat(end);
                        bmaskA = lonmask&latmask;
                        plot(LON(bmaskA), LAT(bmaskA),'.y');
                       
                    end
                    plot(cluster_pnts0(rm_ind,1), cluster_pnts0(rm_ind,2),'xc','markersize',8);
                    

                   
                end
                
                if ~isempty(ind)
                    num_indv_wblob0 = num_indv_wblob;
                    % mask the curent maxPntIdx out:
                    BboxSize(rm_ind) = NaN;
                    cluster_pnts(rm_ind,:) = NaN;
                else
                    search_flag=false;                                 % stop the search loop
                    disp(['final number of blob = ', num2str(num_indv_wblob)]);
                    if ~isempty(ind)
                        indv_blobIndx = [indv_blobIndx, PntIdx(ind')];
                    end
                end
                
                
                    
               
            end
            
            if debug
                % plot them out to confirm:
                figure(1); hold on
                subplot(2,1,2);
                
                plot(blob_center(indv_blobIndx,1),blob_center(indv_blobIndx,2),'pm','markersize',12,'markerFaceColor','m');
                hold on;
                patch(bbox_vertices_x(:,indv_blobIndx), bbox_vertices_y(:,indv_blobIndx),'r', 'FaceColor','None');
               
            end

            blobIndices = [blobIndices, indv_blobIndx];
            
        end                 % end searching for each cluster.
        
        Nblobs = length(blobIndices);
        % use index in blobIndices to select the blob statistics out:
        blobs_sel.GeoLocs = blob_center(blobIndices,:);
        blobs_sel.BoundingBoxSize =  blob_BBoxSize(blobIndices,:);
        % save blob mask:
        blob_Images =  [blobs.blob_image];
        
        for i = 1:length(blobIndices)
            blobs_sel.blobImage{i} = blob_Images{blobIndices(i)};
            blobs_sel.blobPixelIds{i} = blob_PixelIds_all{blobIndices(i)};
        end
        % sabe blobImage coord:
        blobImageCoord = [blobs.blob_image_GeoCoord];
        blobs_sel.blobImageCoord = blobImageCoord(blobIndices); 

        stats_tmp = [blobs.stats_selected];
        stats_vars = fieldnames(blobs(1).stats_selected);
        for i = 1:length(stats_vars)
            SN = stats_vars{i};
            %stats_tmp = [];
            % gather up statistics associated with all the blobs:
            % only take the variables that has 1 dimention
            if numel(blobs(1).stats_selected.(SN)) == length(blobs(1).stats_selected.(SN)) 
                if ~strcmp(SN, 'maxFeretCoord')
%                     for j = 1:nthres
%                         stats_tmp = [stats_tmp, blobs(j).stats_selected.(SN)];     % need to see if vertcat is needed.
%                     end
                      var_tmp = [stats_tmp.(SN)];
                      if ~iscell(var_tmp)
                          blobs_sel.stats_selected.(SN) = var_tmp(blobIndices);
                      else
                          for ii = 1:length(blobIndices)
                              blobs_sel.stats_selected.(SN) = var_tmp{blobIndices(ii)};
                          end
                      end
                end
            end
        end
        
 
return