function cldfrac = compute_spatial_cloud_fraction_in_coarse_grid_parallel(ngrid, ncldmask, cgrid)
% Purpose: compute cloud fraction to describe cloudiness at a coarser grid;
% Inputs:
%  - ngrid (structure) .LON, .LAT (native finer grid that has cloud mask
%  info.)s
%  - ncldmask (array)
%  - cgrid (structure) .LON, .LAT (coraser grid)
% Output: cldfrac (array with (NLat, NLon, Ntime) dimension.
% 

% 1. build pixels centered at grid points:
cgrid.LONm =  0.5*(cgrid.LON(1:end-1,1:end-1) + cgrid.LON(1:end-1,2:end));
cgrid.LATm =  0.5*(cgrid.LAT(1:end-1,1:end-1) + cgrid.LAT(2:end,1:end-1));

ngrid.LONm =  0.5*(ngrid.LON(1:end-1,1:end-1) + ngrid.LON(1:end-1,2:end));
ngrid.LATm =  0.5*(ngrid.LAT(1:end-1,1:end-1) + ngrid.LAT(2:end,1:end-1));
ngrid_LON = ngrid.LON(:); ngrid_LAT = ngrid.LAT(:);

% 2. use ploygin to search for native grid points that fell into the coarse
% pixel
npixel = numel(cgrid.LON) - 2*size(cgrid.LON,2) - 2*(size(cgrid.LON,1)-2);
LONsub = cgrid.LON(2:end-1, 2:end-1);
LATsub = cgrid.LAT(2:end-1, 2:end-1);
LONsubvec = reshape(cgrid.LON(2:end-1, 2:end-1),[],1);
LATsubvec = reshape(cgrid.LAT(2:end-1, 2:end-1),[],1);


res_c = mean(diff(cgrid.LON,1,2),'all'); 
res_n = mean(diff(ngrid.LON,1,2),'all');
c2nratio = res_c^2/res_n^2;

NT = size(ncldmask,3);
cldfrac_vec = NaN(npixel, NT);

tic
hw = waitbar(0, 'converting to cloud frac..');
for it = 1:NT
    ncldmask_tmp = reshape(ncldmask(:,:,it),[],1);
    
    if sum(ncldmask_tmp,'all')==numel(ncldmask_tmp)
        % all clouds
        cldfrac_vec(:, it) = 1;
    else
        
        cldfrac_vec_tmp = zeros(1, npixel);
        
        %extraPnts_lat = zeros(1,5);
        %extraPnts_lon = zeros(1,5);
        %extraPnts_cldmask = zeros(1,5);
        parfor i = 1:npixel
            %i = j;
            %disp(['i=', num2str(i)])
            ccenx = LONsubvec(i);
            cceny = LATsubvec(i);
            
            % get the vertice for this pixel:
            vx = [ccenx-0.5*res_c, ccenx-0.5*res_c, ccenx+0.5*res_c, ccenx+0.5*res_c, ccenx-0.5*res_c];
            vy = [cceny-0.5*res_c, cceny+0.5*res_c, cceny+0.5*res_c, cceny-0.5*res_c, cceny-0.5*res_c];
            
            srad = 1.8;
            r = (ngrid_LON - ccenx).^2 + (ngrid_LAT-cceny).^2;
            rmask = r<=srad*res_c;
            
            ncldmask_nearby = ncldmask_tmp(rmask);
            ngrid_LON_nb = ngrid_LON(rmask);
            ngrid_LAT_nb = ngrid_LAT(rmask);
            
            
            % if all cloudy or all clear in the search radius, then no need to compute cloud
            % fraction in the following manner:
            
            if sum(ncldmask_nearby) == numel(ncldmask_nearby)
                % all cloudy:
                cldfrac_vec_tmp(i) = 1.0;
                
            elseif sum(ncldmask_nearby)==0
                cldfrac_vec_tmp(i) = 0.0;
                
            else
            
            [in, on] = inpolygon(double(ngrid_LON_nb), double(ngrid_LAT_nb), double(vx),double(vy));
            all = in|on;
            
            lon_sel = ngrid_LON_nb(all); lat_sel=ngrid_LAT_nb(all);
            ncldmask_sel = ncldmask_nearby(all);
            
            Aout =0;
            %clear extraPnts_cldmask
            if length(lon_sel)<=6
                %disp(['finding extra points = '  num2str(9-length(lon_sel)) ]);
                % need to include extra points to compute cloud coverage;
                % determine which direction to include the extra points:
                dist2top = min(abs(lat_sel - vy(2)));
                dist2bot = min(abs(lat_sel - vy(1)));
                dist2lft = min(abs(lon_sel - vx(1)));
                dist2rgt = min(abs(lon_sel - vx(3)));
                
                direction_matrix = [-1, 1];
                lat_start = [min(lat_sel), max(lat_sel)];
                lon_start = [min(lon_sel), max(lon_sel)];
                
                add_row = ismembertol([dist2bot, dist2top], res_n, 1e-3);
                add_col = ismembertol([dist2lft, dist2rgt], res_n, 1e-3);
                
                added_lon = []; added_lat = [];
                if any(add_row)
                    added_lon = unique(lon_sel);
                    added_lat = repmat(lat_start*add_row'+ direction_matrix*add_row'*res_n, size(added_lon));
                end
                
                added_lon2 = []; added_lat2 = [];
                if any(add_col)
                    added_lat2 = unique([lat_sel; added_lat]);
                    added_lon2 = repmat(lon_start*add_col' + direction_matrix*add_col'*res_n, size(added_lat2));
                    
                end
                %addedPnts_pair = [];
                A = [added_lon; added_lon2];
                B = [added_lat; added_lat2];
                addedPnts_pair = horzcat(A,B);
                
                addedPnts_pair = unique(addedPnts_pair,'row');
                
                extraPnts_lon = []; %zeros(length(addedPnts_pair),1);
                extraPnts_lat = []; %zeros(length(addedPnts_pair),1);
                extraPnts_lon = addedPnts_pair(:,1);
                extraPnts_lat = addedPnts_pair(:,2);
                
                % find cldmask at the extra points:
                %extraPnts_cldmask=zeros(length(addedPnts_pair),1);
                box = [];
                for ip = 1:length(extraPnts_lon)
                    % not longer 2D, and so, need to find another way.. to find
                    % this..
                    %                [~, mid1] = min(abs(extraPnts_lon(ip)-ngrid_LON));
                    %                [~, mid2] = min(abs(extraPnts_lat(ip)-ngrid_LAT));
                    %
                    %                extraPnts_cldmask(ip) = ncldmask_tmp(mid2, mid1);
                    locmask = ismembertol(ngrid_LON_nb, extraPnts_lon(ip),1e-4) & ismembertol(ngrid_LAT_nb,extraPnts_lat(ip),1e-4)
                    box = [box, ncldmask_nearby(locmask)];
                    
                end
                %disp(num2str(box))
                extraPnts_cldmask = box;
                
                %if exist('extraPnts_lon','var')
                corner_grid = ismembertol(extraPnts_lon, vx, 1e-4) & ismembertol(extraPnts_lat, vy, 1e-4);
                if any(corner_grid)
                    Aout = Aout + sum(extraPnts_cldmask(corner_grid).*0.25*1);
                end
                
                
                border_grid = (ismembertol(extraPnts_lon, vx, 1e-4) | ismembertol(extraPnts_lat, vy, 1e-4)) & ~corner_grid;
                if any(border_grid)
                    Aout = Aout + sum(extraPnts_cldmask(border_grid).* 0.5*1);
                end
                
                
            end
            
            if length(lon_sel)==9 || (length(lon_sel) + length(extraPnts_lon))==9
                
                % double check for on points
                %         on_locs = [];
                %         % use ismember to replace this loop:
                %         for m =1:length(lon_sel)
                %             dtmp = min(abs(lon_sel(m) - vx));
                %             dtmp2 = min(abs(lat_sel(m) - vy));
                %             if dtmp<1e-4 || dtmp2<1e-4
                %                 on_locs= [on_locs, m];
                %             end
                %         end
                on_locs = find( (ismembertol(lon_sel, vx,1e-4) | ismembertol(lat_sel, vy, 1e-4)) )
                %disp(['on_locs=', num2str(on_locs)]);
               % disp(['lon_sel=', num2str(lon_sel)]);
%                if isempty(on_locs)
%                    disp('bug!')
%                else
%                    disp(num2str(on_locs'))
%                end
                
                OnPnts_lon = lon_sel(on_locs);
                OnPnts_lat = lat_sel(on_locs);
                OnPnts_cldmask = ncldmask_sel(on_locs);
                
                lon_sel(on_locs) = [];
                lat_sel(on_locs) = [];
                ncldmask_sel(on_locs)=[];
                
                %disp(['ncldmask_sel=', num2str(ncldmask_sel)]);
                
                InPnts_lon = lon_sel;
                InPnts_lat = lat_sel;
                InPnts_cldmask = ncldmask_sel;
                
                %                 figure(11)
                %                 plot(InPnts_lon, InPnts_lat,'*k');
                %                 hold on
                %                 plot(OnPnts_lon, OnPnts_lat,'*k');
                %                 plot(extraPnts_lon, extraPnts_lat,'*r');
                %                 plot(ccenx, cceny,'or');
                %                 plot(vx,vy,'-r');
                %                 pause
                
                
                Ain = sum(InPnts_cldmask.* 1);    % area of native grid pixels that are completely within the coarser pixel, weighted by cloud mask;
                %disp(['Ain=', num2str(Ain)]);
                
                if ~isempty(on_locs)
                    corner_grid = ismembertol(OnPnts_lon, vx, 1e-4) & ismembertol(OnPnts_lat, vy, 1e-4);
                    if any(corner_grid)
                        Aout = Aout + sum(OnPnts_cldmask(corner_grid).*0.25*1);
                    end
                    
                    
                    border_grid = (ismembertol(OnPnts_lon, vx, 1e-4) | ismembertol(OnPnts_lat, vy, 1e-4)) & ~corner_grid;
                    if any(border_grid)
                        Aout = Aout + sum(OnPnts_cldmask(border_grid).* 0.5*1);
                    end
                end
                
                
                
%                 disp(['Ain=' num2str(Ain)]);
%                 disp(['Aout=' num2str(Aout)]);
                %pause(0.2)
                cldfrac_vec_tmp(i) = double((Ain + Aout))/c2nratio;
                
                %         if cldfrac_vec(i,it)-1>1e-2
                %             length(InPnts.lon)
                %             length(OnPnts.lon)
                %             pause
                %         end
            else
                disp(['it = ' num2str(it) '; ipixel = ' num2str(i) '; total native grid point /=9!']);
                cldfrac_vec_tmp(i) = NaN;
            end
            end
            
            
        end
        
        cldfrac_vec(:,it) = cldfrac_vec_tmp;
        
        
    end
    
    if mod(it, floor(NT/10))<1e-2
        waitbar(it/NT, hw, [num2str(it/NT*100) '% completed..']);
    end
    
    
    
end
cldfrac = reshape(cldfrac_vec, size(LONsub,1), size(LONsub,2),NT);
toc

return
