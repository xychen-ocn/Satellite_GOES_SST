function windInfo = compute_winddiv_over_features(CCMPv2_dataFN, t, search_box, blobsIn, cutoff_thres, varargin)
% Purpose: this script is used to compute wind divergence over features
% from input wind data (CCMPv2 or ERA5)
switch nargin 
    case 5
        masktype = 'circle';
    case 6
        masktype = 'square';
end
  

if ~exist(CCMPv2_dataFN)
    disp('no such file; check input data filename again!')
    return

else
    % read in the daily averaged wind data: (0.25 degree, different
    % resolution from the SST (5km);
    wind_data = read_netCDF_into_matlab_structure(CCMPv2_dataFN,'DimInOrder',{'longitude','latitude','time'});
    
    % compute wind divergence using central difference method:
    % ---- compute wind divergence: ---- %
    [LON_wnd, LAT_wnd] = meshgrid(wind_data.longitude, wind_data.latitude);
    
    dx = diff(LON_wnd,1, 2);
    dx = [zeros(size(dx,1),1)  dx];
    
    dy = diff(LAT_wnd,1,1);
    dy = [zeros(1,size(dy,2)); dy];
    
    XX_L4 = 111*cosd(LAT_wnd).*cumsum(dx,2);     % km
    YY_L4 = 111.*cumsum(dy,1);                        % km
    
    % compute the divergence for this thing:
    tidx = find(floor(wind_data.time_num) == floor(t));
    uwnd = wind_data.uwnd(:,:,tidx)';    % need to permute to have lat in the first dim;
    vwnd = wind_data.vwnd(:,:,tidx)';
    winddiv = divergence(XX_L4, YY_L4, uwnd, vwnd);
    
    % add spatial filtering:
    winddiv_LS =estimate_LSG_SST(LON_wnd(1,:), LAT_wnd(:,1), winddiv, 'method','spectrum','CutoffScale', cutoff_thres, 'checkflag',false);    % the size is different here.
    winddiv_HF = winddiv - winddiv_LS;
    
    % now extract the wind divergence cutout for every individual feature:
    BndBox.cen = blobsIn.GeoLocs;
    BndBox.Size = blobsIn.BoundingBoxSize;
    EqvDiam_dgr = blobsIn.stats_selected.EqvDiam .* blobsIn.L4_xres/111E3;  % in degrees.
    search_RadRatio = search_box.sRadiusRatio;
    cutout_extend = search_RadRatio;
    
    nblobs = length(BndBox.cen);
    
    for ib = 1:nblobs
        xb = BndBox.cen(ib,1); yb = BndBox.cen(ib,2);
        
        if strcmp(masktype, 'circle') % do the following
            
            L = max(BndBox.Size(ib,:))*search_RadRatio;
            
            dist = sqrt((LON_wnd - xb).^2 + (LAT_wnd-yb).^2);
            smask = dist<=L*1.5;
            search_lon = LON_wnd(smask); search_lat = LAT_wnd(smask);
            
            search_lon = LON_wnd(smask);
            search_lat = LAT_wnd(smask);
            
            % extract the pixels from winddiv:
            windInfo(ib).winddiv= winddiv(smask);
            windInfo(ib).winddiv_highfreq = winddiv_HF(smask);
           
            % extract uwnd and vwind separately.
            uwnd_cutouts2D = uwnd(smask);
            vwnd_cutouts2D = vwnd(smask);
            
            
            
        elseif strcmp(masktype, 'square')
            
            
            % 1. extract
            L = EqvDiam_dgr(ib)/2*search_RadRatio;
         
            
            %% under construction:
            % use the same way to cut out cloudiness data as the SST data for
            % comparison.
            LX = BndBox.Size(ib,1);
            LY = BndBox.Size(ib,2);
            
            ExtRad = EqvDiam_dgr(ib)/2*cutout_extend ;
            
            lon0 = xb - LX/2; %- cutout_extend;   % index of the end of the box.
            lonN = xb + LX/2; %+ cutout_extend;
            
            lat0 = yb - LY/2; %- cutout_extend;
            latN = yb + LY/2; %+ cutout_extend;
            
            % create mask:
            cutout_mask = false(size(winddiv));
            cutout_mask((LON_wnd>(lon0-ExtRad))&(LON_wnd<(lonN+ExtRad)))= true;
            cutout_mask((LAT_wnd<(lat0-ExtRad))|(LAT_wnd>(latN+ExtRad)))= false;
            
            latslice = LAT_wnd(:,1);
            nrow = length(find((latslice>=lat0-ExtRad)&(latslice<=latN+ExtRad)));
            
            % extract the pixels from winddiv:
            wndiv_cutouts = winddiv(cutout_mask);
            wndiv_cutouts2D = reshape(wndiv_cutouts, nrow, []);
            windInfo(ib).winddiv = wndiv_cutouts2D;
            
            % extract uwnd and vwind separately.
            uwnd_cutouts = uwnd(cutout_mask);
            vwnd_cutouts = vwnd(cutout_mask);
            uwnd_cutouts2D = reshape(uwnd_cutouts, nrow, []);
            vwnd_cutouts2D = reshape(vwnd_cutouts, nrow, []);
            
            % get lon and lat:
            search_lonvec = LON_wnd(cutout_mask); search_lon = reshape(search_lonvec, nrow, []);
            search_latvec = LAT_wnd(cutout_mask); search_lat = reshape(search_latvec, nrow, []);
            
        end
        
        
        windInfo(ib).time = t;
        windInfo(ib).uwnd = uwnd_cutouts2D;
        windInfo(ib).vwnd = vwnd_cutouts2D;
        windInfo(ib).lon = search_lon;
        windInfo(ib).lat = search_lat;
                
    end
    
    % save data in certain format for mapping to a new coordinate.
    
    
end

return