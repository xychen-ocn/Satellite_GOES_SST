function cldfrac = compute_spatial_cloud_fraction_in_coarse_grid(ngrid, ncldmask, cgrid)
% Purpose: compute cloud fraction to describe cloudiness at a coarser grid;
% Inputs:
%  - ngrid (structure) .LON, .LAT (native finer grid that has cloud mask
%  info.)
%  - ncldmask (array)
%  - cgrid (structure) .LON, .LAT (coraser grid)
% Output: cldfrac (array with (NLat, NLon, Ntime) dimension.
% 

% 1. build pixels centered at grid points:
cgrid.LONm =  0.5*(cgrid.LON(1:end-1,1:end-1) + cgrid.LON(1:end-1,2:end));
cgrid.LATm =  0.5*(cgrid.LAT(1:end-1,1:end-1) + cgrid.LAT(2:end,1:end-1));

ngrid.LONm =  0.5*(ngrid.LON(1:end-1,1:end-1) + ngrid.LON(1:end-1,2:end));
ngrid.LATm =  0.5*(ngrid.LAT(1:end-1,1:end-1) + ngrid.LAT(2:end,1:end-1));

% 2. use ploygin to search for native grid points that fell into the coarse
% pixel
npixel = numel(cgrid.LON) - 2*size(cgrid.LON,2) - 2*(size(cgrid.LON,1)-2);
cgrid.LONsub = cgrid.LON(2:end-1, 2:end-1);
cgrid.LATsub = cgrid.LAT(2:end-1, 2:end-1);

res_c = mean(diff(cgrid.LON,1,2),'all'); 
res_n = mean(diff(ngrid.LON,1,2),'all'); 

NT = size(ncldmask,3);
cldfrac_vec = NaN(npixel, NT);

tic
hw = waitbar(0, 'converting to cloud frac..');
for it = 1:NT
    ncldmask_tmp = ncldmask(:,:,it);
    if sum(ncldmask_tmp,'all')==numel(ncldmask_tmp)
        % all clouds
        cldfrac_vec(:, it) = 1;
    end
    
   
    parfor i = 1:npixel
        
        ccenx = cgrid.LONsub(i); cceny = cgrid.LATsub(i);
        
        % get the vertice for this pixel:
        vx = [ccenx-0.5*res_c, ccenx-0.5*res_c, ccenx+0.5*res_c, ccenx+0.5*res_c, ccenx-0.5*res_c];
        vy = [cceny-0.5*res_c, cceny+0.5*res_c, cceny+0.5*res_c, cceny-0.5*res_c, cceny-0.5*res_c];
        
        srad = 1.8;
        r = (ngrid.LON - ccenx).^2 + (ngrid.LAT-cceny).^2;
        rmask = r<=srad*res_c;
        ncldmask_nearby = ncldmask_tmp(rmask);
        
        % if all cloudy or all clear in the search radius, then no need to compute cloud
        % fraction in the following manner:
        
        if sum(ncldmask_nearby) == numel(ncldmask_nearby)
            % all cloudy:
            cldfrac_vec(i, it) = 1;
            
        elseif sum(ncldmask_nearby)==0
            cldfrac_vec(i,it) = 0;
            
        else  % start doing the following more complicated calculation:

        
%         lon_mask = ngrid.LON - vx(1)>=0 & ngrid.LON-vx(3)<=0;
%         lat_mask = ngrid.LAT - vy(1)>=0 & ngrid.LAT-vy(2)<=0;
%         all = lon_mask & lat_mask;        
        
        [in, on] = inpolygon(double(ngrid.LON), double(ngrid.LAT), double(vx),double(vy));
        all = in|on;
        
        lon_sel = ngrid.LON(all); lat_sel=ngrid.LAT(all);
        ncldmask_sel = ncldmask_tmp(all);

        clear extraPnts
        if length(lon_sel)<=6
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
           addedPnts_pair = [];
           addedPnts_pair(:,1) = [added_lon; added_lon2];
           addedPnts_pair(:,2) = [added_lat; added_lat2];
           
           addedPnts_pair = unique(addedPnts_pair,'row');
           
           extraPnts.lon = addedPnts_pair(:,1);
           extraPnts.lat = addedPnts_pair(:,2);
           
           % find cldmask at the extra points:
           extraPnts.cldmask=[];
           for ip = 1:length(extraPnts.lon)
               [dtmp, mid1] = min(abs(extraPnts.lon(ip)-ngrid.LON(1,:)));
               [dtmp, mid2] = min(abs(extraPnts.lat(ip)-ngrid.LAT(:,1)));

               extraPnts.cldmask(ip) = ncldmask_tmp(mid2, mid1);
           end
        end

       if length(lon_sel)==9 | (length(lon_sel) + length(extraPnts.lon))==9

        % double check for on points
        on_locs = [];
        for m =1:length(lon_sel)
            dtmp = min(abs(lon_sel(m) - vx));
            dtmp2 = min(abs(lat_sel(m) - vy));
            if dtmp<1e-4 || dtmp2<1e-4
                on_locs= [on_locs, m];
            end
        end
        OnPnts.lon = lon_sel(on_locs);
        OnPnts.lat = lat_sel(on_locs);
        OnPnts.cldmask = ncldmask_sel(on_locs);
        
        lon_sel(on_locs) = [];
        lat_sel(on_locs) = [];
        ncldmask_sel(on_locs)=[];
        InPnts.lon = lon_sel;
        InPnts.lat = lat_sel;
        InPnts.cldmask = ncldmask_sel;
        
if 1==0
        % add visualization to check:
        figure(2);clf
        hold on
        plot(ccenx, cceny, 'or');
        hold on;
        plot(cgrid.LONm, cgrid.LATm,'-r','linewidth',1.25); plot(cgrid.LONm', cgrid.LATm', '-r','linewidth',1.25);
        plot(ngrid.LON(all), ngrid.LAT(all),'*k');
        plot(extraPnts.lon, extraPnts.lat, '*c');

        plot(ngrid.LONm, ngrid.LATm,'--k'); plot(ngrid.LONm', ngrid.LATm','--k');
        axis('square')
        axis([ccenx-5*res_c, ccenx+5*res_c, cceny-5*res_c, cceny+5*res_c]);
        %pause(0.1)
end
        
        Ain = sum(InPnts.cldmask.* res_n.^2);    % area of native grid pixels that are completely within the coarser pixel, weighted by cloud mask;
        
        Aout = 0;

        
        if ~isempty(on_locs)
            corner_grid = ismembertol(OnPnts.lon, vx, 1e-4) & ismembertol(OnPnts.lat, vy, 1e-4);
            if any(corner_grid)
                Aout = Aout + sum(OnPnts.cldmask(corner_grid).*0.25*res_n^2);
            end
            
            
            border_grid = (ismembertol(OnPnts.lon, vx, 1e-4) | ismembertol(OnPnts.lat, vy, 1e-4)) & ~corner_grid;
            if any(border_grid)
                Aout = Aout + sum(OnPnts.cldmask(border_grid).* 0.5*res_n^2);
            end
        end
        
        if exist('extraPnts','var') 
            corner_grid = ismembertol(extraPnts.lon, vx, 1e-4) & ismembertol(extraPnts.lat, vy, 1e-4);
            if any(corner_grid)
                Aout = Aout + sum(extraPnts.cldmask(corner_grid).*0.25*res_n^2);
            end
            
            
            border_grid = (ismembertol(extraPnts.lon, vx, 1e-4) | ismembertol(extraPnts.lat, vy, 1e-4)) & ~corner_grid;
            if any(border_grid)
                Aout = Aout + sum(extraPnts.cldmask(border_grid).* 0.5*res_n^2);
            end
            
        end

        
        cldfrac_vec(i,it) = (Ain + Aout)/res_c^2;
        
%         if cldfrac_vec(i,it)-1>1e-2
%             length(InPnts.lon)
%             length(OnPnts.lon)
%             pause
%         end
       else
           disp(['it = ' num2str(it) 'ipixel = ' num2str(i) '; total native grid point /=9!']);
           cldfrac_vec(i,it) = NaN;
       end
        
        end
    end
    
    
    if mod(it, floor(NT/10))<1e-2
        waitbar(it/NT, hw, [num2str(it/NT*100) '% completed..']);
    end
    
end
cldfrac = reshape(cldfrac_vec, size(cgrid.LONsub,1), size(cgrid.LONsub,2),NT);
toc
return
