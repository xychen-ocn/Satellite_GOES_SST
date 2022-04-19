function [ave_windInfo, search_box] = get_ERA5_windInfo_over_features(t, blobsIn, search_RadRatio, ERA5dataFN, varargin)
% purpose: get area-averaged ERA5 10-m surface wind; Area is determined by
% Bounding box and the search radius;
% Input: 
%          t : datenumber, to extract ERA5 wind at specified time. 
%      BndBox: bounding box information (center location, width, height)
%      search_RadRatio: x times the width and height of the bounding box
%      ERA5dataFN: netCDF file that stores ERA5 surface wind.

% read additional inputs if provided:
switch nargin
    case 4
        zwind = 10;
    case 5 
        zwind = varargin{1};
end



ERA5_datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST/era5_data';
ERA5_absFN = [ERA5_datadir filesep ERA5dataFN];

if ~exist(ERA5_absFN)
    disp('no such ERA5 file');
    return
else
    ERA5 = read_netCDF_into_matlab_structure(ERA5_absFN);         % worked fine.
    [LON, LAT] = meshgrid(ERA5.lon-360, ERA5.lat);
        
    [~,tid] = min(abs(ERA5.time_num -t));
    disp(datestr(ERA5.time_num(tid)));
    
    if zwind==10
        uwnd = ERA5.uwnd_10m(:,:,tid)';
        vwnd = ERA5.vwnd_10m(:,:,tid)';
        
    else  % interpolate wind at the request level zwind
        % sudo code:
        % ERA5.uwnd dims: lon, lat, lev, time
        % I need to convert geopotential height to height first;
        %zvec = 1:size(ERA5.hgt,3);
       % [LON3D, LAT3D, Z3D] = meshgrid(ERA5.lon-360, ERA5.lat, zvec);
       % [xq, yq, zq] = meshgrid(ERA5.lon-360, ERA5.lat, zwind:zwind+1);
        
        hgt = permute(ERA5.hgt(:,:,:,tid), [2,1,3]);
        ERA5_u = permute(ERA5.uwnd(:,:,:, tid), [2,1,3]);
        ERA5_v = permute(ERA5.vwnd(:,:,:, tid), [2,1,3]);
        
        uwnd = zeros(size(ERA5_u,1), size(ERA5_u,2));
        vwnd = uwnd;
        for j = 1:length(ERA5.lat)
            for i = 1:length(ERA5.lon)
                uwnd(j,i) = interp1(squeeze(hgt(j,i,:)), squeeze(ERA5_u(j,i,:)), zwind);
                vwnd(j,i) = interp1(squeeze(hgt(j,i,:)), squeeze(ERA5_v(j,i,:)), zwind);
            end
        end

        
    end
    
    
    
    % compute divergence:
    wnd_div_all = compute_div(ERA5.lon, ERA5.lat,uwnd, vwnd);
    
    [BVX, BVY] = get_boundingbox_vertices(blobsIn);
    BndBox.cen = blobsIn.GeoLocs;
    BndBox.Size = blobsIn.BoundingBoxSize;              % LX (width); LY (height);
    
    % get data point within and on the search area: (square with the size
    % the same length as the longest dimension of the bounding box) 
    %% note: in Park et al. (2006) the wind vector cells are either 11x11 or 13x13, from a roughly squared region.)
    for ib = 1:size(BndBox.cen,1)
        BlobMask = blobsIn.blob_image{ib};
        BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);

        
        L = max(BndBox.Size(ib,:))*search_RadRatio;
        
        xb = BndBox.cen(ib,1); yb = BndBox.cen(ib,2);
        
        if 1==0
        % define the search sqaure:
        square_vx = [xb-0.5*L, xb-0.5*L, xb+0.5*L, xb+0.5*L, xb-0.5*L];
        square_vy = [yb-0.5*L, yb+0.5*L, yb+0.5*L, yb-0.5*L, yb-0.5*L];
        
        [in, on] = inpolygon(LON, LAT, square_vx, square_vy);
         
        % find mean wind directionin the search area;
        mask = in | on;
        end
        
        
        
        % define with radii instead:
        dist = sqrt((LON - xb).^2 + (LAT-yb).^2);
        mask = dist<=L;
        
        if isempty(find(mask==1))
            
            disp('no ERA5 grid at the search location');
            % make local grid at ERA5 resolution for interpolation:
            gridx = [min(square_vx): 0.25: max(square_vx)];
            gridy = [min(square_vy): 0.25: max(square_vy)];  
            
            [GX, GY] = meshgrid(gridx, gridy);
            
            %uwnd_vertices = interp2(ERA5.lon-360, ERA5.lat, uwnd, [square_vx(1:4), xb], [square_vy(1:4), yb]);
            %vwnd_vertices = interp2(ERA5.lon-360, ERA5.lat, vwnd, [square_vx(1:4), xb], [square_vy(1:4), yb]);
            uwnd_intp = interp2(ERA5.lon-360, ERA5.lat, uwnd, GX, GY);
            vwnd_intp = interp2(ERA5.lon-360, ERA5.lat, vwnd, GX, GY);
            mean_u = mean(uwnd_intp,'omitnan');
            mean_v = mean(vwnd_intp,'omitnan');
            
        else
            
            mean_u = mean(uwnd(mask),'omitnan');
            mean_v = mean(vwnd(mask), 'omitnan');
            
        end
        
        mean_wnddir = atan2(mean_v, mean_u)*180/pi;
        mean_wndspd = sqrt(mean_u.^2 + mean_v.^2);
       
        
       % disp(['mean_wnddir 00:' num2str(mean_wnddir)])
        % define bounding box in a coordinate where the mean wind is aligned
        % with the x-axis. (positive x-direction: downwind)
        th = atan(mean_v/mean_u);  % [-pi/2, pi/2];
        LX_dw = BndBox.Size(ib,1)*cos(th);
        LY_dw = BndBox.Size(ib,2)/cos(th);
        
        BBox_newCoord_vx = [-0.5*LX_dw, 0.5*LX_dw, 0.5*LX_dw, -0.5*LX_dw, -0.5*LX_dw];
        BBox_newCoord_vy = [-0.5*LY_dw, -0.5*LY_dw, 0.5*LY_dw, 0.5*LY_dw, -0.5*LY_dw];
        BBox_newCoord = vertcat(BBox_newCoord_vx, BBox_newCoord_vy);
        
        R = rotz(mean_wnddir);
        
        BBoxRot_in_cartCoord = R*BBox_newCoord + BndBox.cen(ib,:)';
        
        
        %% visualize two bounding box to check:
        if 1==0
        figure(10); hold on;
        % plot blob element
        pcolor(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask)); shading flat;
        hold on;
        plot(BVX(:,ib), BVY(:,ib),'-','marker','o','color','b');
        hold on;
        
        plot(BBoxRot_in_cartCoord(1,:), BBoxRot_in_cartCoord(2,:), '-r','marker','*');
       % pause
        end
        
        BndBox_dw.Size= [LX_dw, LY_dw];
        BndBox_dw.cen = BndBox.cen;
        BndBox_dw.vertices = BBoxRot_in_cartCoord;
        
        search_box.BndBox_dw(ib) = BndBox_dw;
        
       
        %% if looked okay, then extract wind_divergence for the search region:
        % how to store the information is a problem. do it first and worry
        % about it later;
        if 1==0
        % define the search square again (in the new downwind coordinate)
        SBox_newCoord_vx = [-0.5*LX_dw, 0.5*LX_dw, 0.5*LX_dw, -0.5*LX_dw, -0.5*LX_dw].*search_RadRatio;
        SBox_newCoord_vy = [-0.5*LY_dw, -0.5*LY_dw, 0.5*LY_dw, 0.5*LY_dw, -0.5*LY_dw].*search_RadRatio;
        SBox_newCoord = vertcat(SBox_newCoord_vx, SBox_newCoord_vy);
        
        R = rotz(-mean_wnddir);
        
        SBoxRot_in_cartCoord = inv(R)*SBox_newCoord + BndBox.cen(ib,:)';
        %search_box.SBox_dw_vertices{ib} = SBoxRot_in_cartCoord;
        
        [in, on] = inpolygon(LON, LAT, SBoxRot_in_cartCoord(1,:), SBoxRot_in_cartCoord(2,:));
        
        % extract wind divergence in the search area;
        mask = in | on;
        
        
        mean_u = mean(uwnd(mask),'omitnan');
        mean_v = mean(vwnd(mask), 'omitnan');
        
        mean_wnddir = atan2(mean_v, mean_u)*180/pi;
        mean_wndspd = sqrt(mean_u.^2 + mean_v.^2);
        end
       % disp(['mean_wnddir 01:' num2str(mean_wnddir)])
        
        % update again;
        ave_windInfo(ib).u = mean_u;
        ave_windInfo(ib).v = mean_v;
        ave_windInfo(ib).wnddir = mean_wnddir;
        ave_windInfo(ib).wndspd = mean_wndspd;
        
            
        %% 2x bigger search area to save u,v wind information:
        
        if 1==0
        SBox_newCoord_vx = [-0.5*LX_dw, 0.5*LX_dw, 0.5*LX_dw, -0.5*LX_dw, -0.5*LX_dw].*search_RadRatio*1.5;
        SBox_newCoord_vy = [-0.5*LY_dw, -0.5*LY_dw, 0.5*LY_dw, 0.5*LY_dw, -0.5*LY_dw].*search_RadRatio*1.5;
        SBox_newCoord = vertcat(SBox_newCoord_vx, SBox_newCoord_vy);
        
        R = rotz(mean_wnddir);
        
        SBoxRot_in_cartCoord = R*SBox_newCoord + BndBox.cen(ib,:)';
        search_box.SBox_dw_vertices{ib} = SBoxRot_in_cartCoord;
        
        [in, on] = inpolygon(LON, LAT, SBoxRot_in_cartCoord(1,:), SBoxRot_in_cartCoord(2,:));
        
        % extract wind divergence in the search area;
        mask = in | on;
        end
        
        dist = sqrt((LON - xb).^2 + (LAT-yb).^2);
        mask = dist<=L*1.5;
        
        search_box.lon{ib} = LON(mask);
        search_box.lat{ib} = LAT(mask);
        search_box.wspd{ib} = sqrt(uwnd(mask).^2 + vwnd(mask).^2);
        search_box.wdir{ib} = atan2(vwnd(mask), uwnd(mask))*180/pi;
        
          
        
       % end
        
        
    end
            search_box.sRadiusRatio = search_RadRatio;

end

return