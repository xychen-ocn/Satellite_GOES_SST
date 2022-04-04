function cloudiness = map_cloudiness_to_newCoord(ave_windInfo, cloudiness, blobsIn)
% purpose: map the data in the cartesian coordinate to a feature centric,
%          trade-wind aligned, and normalized coordinate.
%          how to normalize the coordinate? (same as Park et al. 2006);
%  three coordinate: 
%     i. cartesian (x,y) -> all input coord by default;
%    ii. ellipse-centric normalized coordinate (xE, yE);
%   iii. wind-aligned coordinate (xA, yA);
%
% transform coordinates of data points extracted in the search region
% (search_box.lon, search_box.lat; cloudiness.lon, cloudiness.lat);

% Inputs:

% Outputs:

xdir = [ave_windInfo.wnddir];
xdir(xdir<0)= xdir(xdir<0)+360;

nblobs = length(blobsIn.GeoLocs);

for ib = 1:nblobs
    BlobMask = blobsIn.blob_image{ib};
    BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);
    
    bx = blobsIn.GeoLocs(ib,1); by = blobsIn.GeoLocs(ib,2);
    cldm.lon = cloudiness(ib).lon; cldm.lat = cloudiness(ib).lat;
    
    % 1. establish feature centric cartesian coordinate
    x_cc = (cldm.lon - bx)*111E3.*cosd(cldm.lat); y_cc = (cldm.lat - by)*111E3;
    if ~isrow(x_cc)
        x_cc = x_cc'; y_cc = y_cc';
    end
    
    %% to-be updated: 
    % add option here to use circle-centric, normalized coordinate. (try
    % using function)
    % 2. establish ellipse-centric, normalized coordinate:
    theta_ori = blobsIn.stats_selected.maxFeretAng(ib);
    if theta_ori<0
        theta_ori = theta_ori+180;
    end
    majA = blobsIn.stats_selected.MajAxisLen(ib)/2 * blobsIn.L4_xres;               % in pixel, --> convert to km 
    minB = blobsIn.stats_selected.MinAxisLen(ib)/2 * blobsIn.L4_xres;
    
    % direction of xE is aligned with the major axis direction
    % (orientation);
    R_E = rotz(-theta_ori);
    
    EllipseCoord = R_E*[x_cc; y_cc];        % rotated coordinate 
    
    xE = EllipseCoord(1,:);  yE = EllipseCoord(2,:);
    
    % normalization by the ellipse boundary (equiv. boundary of the warm features)
    Rp = sqrt(xE.^2 + yE.^2);
    tan_thp = yE./xE;      
    
    xbp = (majA*minB)./sqrt(minB^2 + majA^2.*tan_thp.^2); 
    ybp = xbp.*tan_thp;
    
    Rbp = sqrt(xbp.^2 + ybp.^2);
    
    Rn = Rp./Rbp;
    thp = atan2(yE, xE);
    xEn = Rn.*cos(thp);
    yEn = Rn.*sin(thp);
       
    
    % 3. rotate the ellipse-centric, normalized coordinate such that its y axis aligns with
    % the wind direction. 
    theta_wind= xdir(ib);
    theta_WE = theta_wind-theta_ori;

    theta_AE = theta_WE - 90;
    
    R_AE = rotz(-theta_AE);
    
    WindAligned_Coord = R_AE * [xEn;yEn];
    
    cloudiness(ib).WindAligned_Coord = WindAligned_Coord;
    cloudiness(ib).mean_wspd = ave_windInfo.wndspd;
    
    % 4. visualize results to double check:
    if ~isnan(cloudiness(ib).cloudfreq)
        if 1==0
        figure(12)
        scatter(cldm.lon, cldm.lat, 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
        hold on
        contour(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask),'-r');
        hl = plot_ellipse(gca, bx, by, majA/111E3, minB/111E3,theta_ori);
        hold off
        %
        
        figure(11); clf;
        subplot(1,3,1);
        %scatter(cldm.lon, cldm.lat, 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
        scatter(x_cc, y_cc, 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
        hold on;
        % plot ellipse:
        hl = plot_ellipse(gca, 0,0, majA, minB,theta_ori);
        % plot wind direction:
        quiver(0,0, ave_windInfo(ib).u*10E3, ave_windInfo(ib).v*10E3,'w', 'linewidth',2);
        axis('square');
        colormap(jet); colorbar;
        caxis([0, 1])
        
        subplot(1,3,2);
        scatter(xE, yE, 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
        hold on
        hl = plot_ellipse(gca, 0,0, majA, minB,0);
        quiver(0,0, cosd(theta_WE)*6*10E3, sind(theta_WE)*6*10E3,'w', 'linewidth',2);
        hold on;
        axis('square');
        colormap(jet); colorbar;
        caxis([0, 1])
        
        % doesn't seem correct??
        subplot(1,3,3);
        scatter(WindAligned_Coord(1,:), WindAligned_Coord(2,:), 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
        hold on;
        circle(0,0,1,'w',2);
        xlabel('fractional distance (cross-wind)');
        ylabel('fractional distacne (down-wind');
        axis('square');
        colormap(jet); colorbar;
        caxis([0, 1])
        end
        %pause
    end
    


    %pause(0.1)
end
return

