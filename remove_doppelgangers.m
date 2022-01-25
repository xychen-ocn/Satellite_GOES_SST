function [blobs_sel, Nblobs]=remove_doppelgangers(blobs)
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
%

blob_cenlon_tmp = []; blob_cenlat_tmp=[];
LXs_tmp = []; LYs_tmp = [];

nthres = length(blobs);

%% 1. aggregate all the blobs caputred by n dT thresolds:
    % a. for given dT_thres i, compare the center of the warm blobs detected.
    % see if they form a cluster cloud. if so, the center of the cluster is
    % defined as the mean location of the warm blobs, at the given threshould.
    blob_center = vertcat(blobs.GeoLocs);                              % this is now a nblob_tot x 2 matrix
    blob_BBoxSize  = vertcat(blobs.BoundingBoxSize);                        % this matrix is also nblob_tot x 2matrix [LXs, LYs]

    blob_cenlon_tmp = blob_center(:,1);
    blob_cenlat_tmp = blob_center(:,2);
    
    LXs_tmp = blob_BBoxSize(:,1);
    LYs_tmp = blob_BBoxSize(:,2);
    
    % b. for given time, the warm blob center can changes too but hopefully not
    % too much, this will allow us to determine the size and shape of the warm
    % blob by using the largest size and shape.
    
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
        axis('equal')
        title('detected blobs with varies \DeltaT thresholds');
        colorbar
        
        % define regions to remove the same blobs
        %cenloc=[blob_cenlon_tmp, blob_cenlat_tmp];  % This is basically blob center
        
        cluster_tag = clusterdata(blob_center, 'Linkage','centroid','MaxClust',10);
        subplot(2,1,2);hold off
        scatter(blob_center(:,1), blob_center(:,2), 20, cluster_tag,'filled');
        colorbar
        colormap(jet);
        hold on;
        patch(bbox_vertices_x, bbox_vertices_y, 'r','FaceColor','none');
        axis('equal')
        title('clustered region to remove warm blobs doppelgangers');
        
        
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
            BoxLx = 0.5*LXs_tmp(Pntmask);
            BoxLy = 0.5*LYs_tmp(Pntmask);
            
            % 2. compute distance between warm blob centroids: check from the
            % centroid with the largest bounding box: compute distance of the rest
            % of the warm blob centroids in the cluster to it.
            disp(['input num. of blobs=' num2str(length(BoxLx))]);
            search_flag=true;
            cluster_pnts = cluster_pnts0;
            num_indv_wblob0 = length(cluster_pnts);
            isearch = 0;
            while search_flag
                isearch=isearch+1;
                [maxL, maxPntIdx] = max(BboxSize);
                deltaX = abs(cluster_pnts(:,1) -cluster_pnts(maxPntIdx,1));
                deltaY = abs(cluster_pnts(:,2) -cluster_pnts(maxPntIdx,2));
                dist = sqrt(deltaX.^2 + deltaY.^2);
                
                maxBoxLx = BoxLx(maxPntIdx);
                maxBoxLy = BoxLy(maxPntIdx);
                
                cond_Lx = deltaX>maxBoxLx; cond_Ly = deltaY>maxBoxLy;
                cond_Ldiag = dist>maxL;
                % the centroid is outside of the biggest (selected) bounding box.
                ind = find(cond_Lx|cond_Ly|cond_Ldiag);     % these are the blobs that are considered different from the maximum blob.
                rm_ind = find(~(cond_Lx|cond_Ly|cond_Ldiag));    % included the current maxPntIdx.
                
                % check if there are other things remo
                
                %rm_ind = find(~cond_Lx& ~cond_Ly)
                figure(1);
                subplot(2,1,2); hold on;
                plot(cluster_pnts0(rm_ind,1), cluster_pnts0(rm_ind,2),'xc','markersize',8);
                % pause
                
                num_indv_wblob = (isearch-1) + 1 + length(ind);
                disp(['num individual blob=', num2str(num_indv_wblob)]);
                indv_blobIndx=[indv_blobIndx, PntIdx(maxPntIdx)];
                
                
                if ~isempty(ind)
                    num_indv_wblob0 = num_indv_wblob;
                    % mask the curent maxPntIdx out:
                    BboxSize(rm_ind) = NaN;
                    cluster_pnts(rm_ind,:) = NaN;
                else
                    search_flag=false;
                    disp(['final number of blob = ', num2str(num_indv_wblob)]);
                    if ~isempty(ind)
                        indv_blobIndx = [indv_blobIndx, PntIdx(ind')];
                    end
                end
                
                
            end
            
            
            % plot them out to confirm:
            figure(1);
            subplot(2,1,2);
            
            plot(blob_center(indv_blobIndx,1),blob_center(indv_blobIndx,2),'pm','markersize',12,'markerFaceColor','m');
            hold on;
            h=patch(bbox_vertices_x(:,blobIndices), bbox_vertices_y(:,blobIndices),'r', 'FaceColor','None');
            

            blobIndices = [blobIndices, indv_blobIndx];
            
        end
        
        Nblobs = length(blobIndices);
        % use index in blobIndices to select the blob statistics out:
        blobs_sel.GeoLocs = blob_center(blobIndices,:);
        blobs_sel.BoundingBoxSize =  blob_BBoxSize(blobIndices,:);
        % save blob mask:
        blob_Images =  [blobs.blob_image];
        for i = 1:length(blobIndices)
            blobs_sel.blobImage{i} = blob_Images{blobIndices(i)};
        end
        % sabe blobImage coord:
        blobImageCoord = [blobs.blob_image_GeoCoord];
        blobs_sel.blobImageCoord = blobImageCoord(blobIndices); 

        
        stats_vars = fieldnames(blobs(1).stats_selected);
        for i = 1:length(stats_vars)
            SN = stats_vars{i};
            stats_tmp = [];
            % gather up statistics associated with all the blobs:
            % only take the variables that has 1 dimention
            if numel(blobs(1).stats_selected.(SN)) == length(blobs(1).stats_selected.(SN)) 
                if ~strcmp(SN, 'maxFeretCoord')
                    for j = 1:nthres
                        stats_tmp = [stats_tmp, blobs(j).stats_selected.(SN)];     % need to see if vertcat is needed.
                    end
                    
                    blobs_sel.stats_selected.(SN) = stats_tmp(blobIndices);
                end
            end
        end
        
 
return