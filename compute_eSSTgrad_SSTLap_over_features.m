function [eSSTgrad, SSTlaplacian] = compute_eSSTgrad_SSTLap_over_features(ave_windInfo, blobsIn)
% Purpose: use ERA5 averaged wind vector to compute the effective SST
% gradient (SST rate of change felt by air parcle)
%
% reference script:
%  explore relation between cloudiness and airTtend.m
%
addpath '/Users/xchen/Documents/GitHub/SAM_LES/matlab'
nblobs = length(ave_windInfo);
for ib = 1:nblobs
    % 1. compute SST gradient
    % should I remove the large scale gradient or not?
    SSTa_tmp = blobsIn.SSTa_cutouts{ib};
    LON = blobsIn.LON_cutouts{ib}; 
    LAT = blobsIn.LAT_cutouts{ib};
    
    lon = LON(1,:)'; lat = LAT(:,1);
    meanlat = mean(lat);
    xvec = [0; cumsum(diff(lon)*111.*cosd(meanlat))];   % km
    yvec = [0; cumsum(diff(lat)*111)];                 % km

    SSTa_grad(ib) = spatial_grad(SSTa_tmp, xvec, yvec);
    SSTlaplacian{ib} = spatial_laplacian(SSTa_tmp, xvec,yvec)*100;     % K/km^2  --> K/100km^2;
    
    % 2. compute effective SST rate of change (multiply by the area-averaged
    % wind vector or the coarse 0.25° wind??
    ERA5_uc = ave_windInfo(ib).u;
    ERA5_vc = ave_windInfo(ib).v;
    eSSTgrad(ib).val = (ERA5_uc.*SSTa_grad(ib).xcomp + ERA5_vc .* SSTa_grad(ib).ycomp)*3.6;  % units: K/h; (1m/s = 3.6km/h)
    eSSTgrad(ib).units = 'K/h';

    % 3. sanity check by plotting out results:  % add contours to make it
    % better? 
if 1==0
    figure(12); clf;
    subplot(1,2,1)
    pcolor(LON, LAT, SSTa_tmp); shading flat;
    hold on
    midlon = mean(lon); midlat =mean(lat);
    quiver(midlon, midlat, ERA5_uc, ERA5_vc,'k');
    colorbar
    colormap(redblue);

    subplot(1,2,2);
    pcolor(LON, LAT, eSSTgrad(ib).val); shading flat;
    hold on
    colorbar
    colormap(redblue); caxis([-4, 4]*10^(-4));
    %pause
end

    
end




return