%function [cloudiness, CF_regridded] = map_cloudiness_to_newCoord(ave_windInfo, cloudiness, blobsIn, varargin)
 function dataIn = map_data_to_newCoord_v2(ave_windInfo,dataIn, blobsIn, varargin)
%  
% purpose: map the data in the cartesian coordinate to a feature centric,
%          trade-wind aligned, and normalized coordinate.
%          how to normalize the coordinate? (same as Park et al. 2006);
%             - two options: a. Park et al. 2006 (ellipse); b. 
%  three coordinate: 
%     i. cartesian (x,y) -> all input coord by default;
%    ii. ellipse/circle-centric normalized coordinate (xE, yE);
%   iii. wind-aligned coordinate (xA, yA);
%
% transform coordinates of data points extracted in the search region
% (search_box.lon, search_box.lat; cloudiness.lon, cloudiness.lat);
%
% Note: Apr 6, 2022 --> generalize this function so that it deals not on ly
% with cloudinss but also with SST and potentially other data;

% Inputs:

% Outputs:

default_xgrd = [-3:0.05:3];
default_ygrd = default_xgrd;
expected_varns={'cloudfreq', 'cloudfrac', 'SST_cutouts', 'SSTa_cutouts', 'eSSTgrad'};

%% parse inputs:
default_option = 'ellipse';
p = inputParser;
addRequired(p,'ave_windInfo',@isstruct);
addRequired(p,'dataIn', @isstruct);
addRequired(p,'blobsIn', @isstruct);
checkString = @(s) any(strcmpi(s, {'ellipse','circle'}));
addParameter(p,'options', default_option, checkString);
addParameter(p,'xgrid', default_xgrd, @isnumeric);
addParameter(p,'ygrid', default_ygrd, @isnumeric);

parse(p, ave_windInfo, dataIn, blobsIn, varargin{:});
ave_windInfo = p.Results.ave_windInfo;
dataIn = p.Results.dataIn;
blobsIn = p.Results.blobsIn;
options = p.Results.options;
xgrd = p.Results.xgrid;       % grid vector for the normalized feature centered, normalized, wind aligned coordinate.
ygrd = p.Results.ygrid;

%% start function:
xdir = [ave_windInfo.wnddir];
xdir(xdir<0)= xdir(xdir<0)+360;

nblobs = length(blobsIn.GeoLocs);

for ib = 1:nblobs
    BlobMask = blobsIn.blob_image{ib};
    BlobMaskCoord = blobsIn.blob_image_GeoCoord(ib);
    
    bx = blobsIn.GeoLocs(ib,1); by = blobsIn.GeoLocs(ib,2);

    % identify  longitude, latitude variables:
    varnames_In = fieldnames(dataIn);
    idx = find(contains(varnames_In,'lon','IgnoreCase',true));
    lon_varn = varnames_In{idx}; 

    idx = find(contains(varnames_In,'lat','IgnoreCase',true));
    lat_varn = varnames_In{idx};

    %if length(dataIn)==nblobs
        data.lon = dataIn(ib).(lon_varn); data.lat = dataIn(ib).(lat_varn);
%     elseif length(dataIn) == 1
%         data.lon = dataIn.(lon_varn){ib}; 
%         data.lat = dataIn.(lat_varn){ib};
%     end
      if isempty(data.lon)
          data.lon = nan(1);
          data.lat = nan(1);
      end
       
    
    % 1. establish feature centric cartesian coordinate
    x_cc = (data.lon - bx)*111E3.*cosd(data.lat); y_cc = (data.lat - by)*111E3;
    
   
        
    
    if numel(x_cc)==length(x_cc) & ~isrow(x_cc)    % a vector
        x_cc = x_cc'; y_cc = y_cc';
    end
    
    if numel(x_cc)~=length(x_cc)    % a matrix
        x_cc = reshape(x_cc,1,[]);    % reshape to be a row vector:
        y_cc = reshape(y_cc,1,[]);
    end
    
   
    %% compute needed dimensions information for each blob;
    D = blobsIn.stats_selected.EqvDiam(ib)* blobsIn.L4_xres;           % units: m
    
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
    
    
    %% to-be updated:
    %% double check if the following code modification is correct..
    % add option here to use circle-centric, normalized coordinate. (try
    % using function)
    % 2. establish ellipse-centric, normalized coordinate:
    if strcmpi(options, 'circle')
        % 1. normalize the coordinate
        Rb = D/2;
        
        % quite similar to the ellipse: first normalize in the polar
        % coordinate, and then transform back to catesian
        %Rp = sqrt(x_cc.^2 + y_cc.^2);
        Rp = sqrt(xE.^2 + yE.^2);
        
        Rn = Rp./Rb;
        %thp = atan2(y_cc, x_cc);
        thp = atan2(yE, xE);
        xEn = Rn.*cos(thp);
        yEn = Rn.*sin(thp);
        
        % 2. rebase the coordinate so that the y-axis is aligned with the
        % wind direction.
        %% update note: Apr 24: instead of make y-axis aligned with wind direction, make it aligned with the SST gradient downwind..
        % this will need to be edited.
        %         theta_wind= xdir(ib);
        %         theta_WE = theta_wind-0;  % delta theta from the x-axis. Here, x-axis is due east.
        %
        %         theta_AE = theta_WE - 90;   %? delta theta from the y-axis.
        %
        %         Rot_AE = rotz(-theta_AE);
        %
        %         WindAligned_Coord = Rot_AE * [xEn;yEn];
        
        
        
        
    elseif strcmpi(options, 'ellipse')
        
        
        
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
        
    end
    % 3. rotate the ellipse-centric, normalized coordinate such that its y axis aligns with
    % the wind direction.
    % rotate to match the downwind SST gradient instead. (north will
    % always corresponds to the maximum effective gradient)
    theta_wind= xdir(ib);
    theta_WE = theta_wind-theta_ori;   % angle between wind and the blob orientation.
    %
    %         theta_AE = theta_WE - 90;
    %
    %         Rot_AE = rotz(-theta_AE);
    %
    %         WindAligned_Coord = Rot_AE * [xEn;yEn];
    
    if abs(theta_WE)>90
        ang_wind_majA = 180-abs(theta_WE);
    else
        ang_wind_majA = abs(theta_WE);
    end
    
    if ang_wind_majA>=45
        % wind vector closer y-axis;
        % no need to rotate y-axis;
        WindAligned_Coord = [xEn; yEn];
    else
        % wind vector closer to x-axis:
        % rotate ellipse major axis to be the y axis:
        
        % -90 or 90 needs to be based on the sign of theta_WE:
        rotdir = sign(cosd(theta_WE));   % Is this a potential problem??
        Rot_AE = rotz(rotdir*90);
        
        WindAligned_Coord = Rot_AE * [xEn; yEn];
    end
        
                
    %---------------- the above are latest modification that needed double check -------- %        
        
%         cloudiness(ib).WindAligned_Coord = WindAligned_Coord;
%         cloudiness(ib).mean_wspd = ave_windInfo.wndspd;
        
  
    
    dataIn(ib).WindAligned_Coord = WindAligned_Coord;
    dataIn(ib).mean_wspd = ave_windInfo.wndspd;
    
    
    % 4. visualize results to double check:
    contain_var = ismember(varnames_In, expected_varns);
    if any(contain_var)
        idx = find(contain_var);
        varn = varnames_In{max(idx)};
        if ~isnan(dataIn(ib).(varn))
            if 1==0
                figure(12)
                if numel(data.lon) == length(data.lon)
                    scatter(data.lon(:), data.lat(:), 50, dataIn(ib).(varn)(:), 'filled','marker','s');
                else
                    pcolor(data.lon, data.lat, dataIn(ib).(varn)); shading flat;
                end
                
                hold on
                contour(BlobMaskCoord.lon, BlobMaskCoord.lat, double(BlobMask),'-r');
                if strcmpi(options, 'ellipse')
                    hl = plot_ellipse(gca, bx, by, majA/111E3, minB/111E3,theta_ori);
                elseif strcmpi(options, 'circle')
                    % plot circle instead:
                    
                end
                hold off
                %pause
                %
                
                figure(11); clf;
                subplot(1,3,1);
                %scatter(cldm.lon, cldm.lat, 50, cloudiness(ib).cloudfreq, 'filled','marker','s');
                scatter(x_cc, y_cc, 50, dataIn(ib).(varn)(:), 'filled','marker','s');
                hold on;
               % if strcmpi(options, 'ellipse')
                    % plot ellipse:
                    hl = plot_ellipse(gca, 0,0, majA, minB,theta_ori);
                %elseif strcmpi(options, 'circle')
                    % plot circle instead:
                    hl = circle(0,0,Rb, 'w',2);
                    
                %end
                % plot wind direction:
                quiver(0,0, ave_windInfo(ib).u*10E3, ave_windInfo(ib).v*10E3,'k', 'linewidth',2);
                axis('square');
                colormap(jet); colorbar;
                title(varn)
                %caxis([0, 1])
                
                
               % if strcmpi(options, 'ellipse')
                    subplot(1,3,2);
                    scatter(xE, yE, 50, dataIn(ib).(varn)(:), 'filled','marker','s');
                    hold on
                    
                    hl = plot_ellipse(gca, 0,0, majA, minB,0);
                    
                    quiver(0,0, cosd(theta_WE)*6*10E3, sind(theta_WE)*6*10E3,'k', 'linewidth',2);
                    hold on;
                    axis('square');
                    colormap(jet); colorbar;
                    %caxis([0, 1])
              %  end
                
                % doesn't seem correct??
                subplot(1,3,3);
                scatter(WindAligned_Coord(1,:), WindAligned_Coord(2,:), 50, dataIn(ib).(varn)(:), 'filled','marker','s');
                hold on;
                circle(0,0,1,'w',2);
                xlabel('fractional distance (cross-wind)');
                ylabel('fractional distacne (down-wind');
                axis('square');
                colormap(jet); colorbar;
                %caxis([0, 1])
                
            end
            %pause
        end
    end

     %%%%%%%%%
     % e2: regrid the output above. (need to be developed again.),
     % CF_regridded will be used directly to explor in different
     % conditions.
     
     % what should be stcIn?  cloudiness is a matlab structure that
     % contains multiple fields with a length equal to the number of
     % features identified from the cloud mask.
     strcIn = dataIn(ib);
     
     % -------------------- To-Follow up ------------------------------- %
     % TO-DO: The following function only applies for cloud masks right now, I
     % need it to work for SST. 
     % need to think about what interpolation method makes more sense.
     dataIn(ib).regridded = grid_data_into_windaligned_coord(strcIn, xgrd, ygrd);    % output size = 1x1; 
     
    
     % check if regrid is successfuly
     if 1==0
     if ~all(isnan(dataIn(ib).regridded.(varn)))%ib ==1 || ib==nblobs
         figure(12);clf;
         pcolor(dataIn(ib).regridded.XX, dataIn(ib).regridded.YY, dataIn(ib).regridded.(varn));
         shading flat;
         colorbar;
         colormap(jet);
         pause(0.2)
     end
     end
     


        
end   % end looping through all features (blobs)
    
    
    
    %pause(0.1)


return

