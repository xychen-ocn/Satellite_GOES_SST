function [SST_LS, SST_maps_local_mean] =  estimate_LSG_SST(lon, lat, SST_maps , varargin)
% input: lon (vector), lat (vector), SST_maps(matrix) (required)
%        optional argument: 'method', 'CutoffScale','checkflag'
default_method = 'planefit';
default_checkflag = true;
default_LScutoff = 1000;    % km
%% build up the parse object to read in element;

p = inputParser;
% define parsing rules:
addRequired(p, 'lon', @isnumeric);
addRequired(p, 'lat', @isnumeric);
addRequired(p,'SST_maps', @isnumeric);
checkString = @(s) any(strcmpi(s, {'planefit','spectrum'}));
addParameter(p,'method', default_method, checkString);
addParameter(p, 'CutoffScale', default_LScutoff, @isnumeric);           % units: km
addParameter(p, 'checkflag', default_checkflag, @islogical);

% ready to parse the inputs:
parse(p, lon, lat, SST_maps,varargin{:});
lon = p.Results.lon;
lat = p.Results.lat;
SST_maps = p.Results.SST_maps;
method = p.Results.method;
LS_cutoff = p.Results.CutoffScale;
checkflag = p.Results.checkflag;


SST_maps_local_mean = mean(SST_maps, 3, 'omitnan');
SST_maps_local_mean0 = SST_maps_local_mean;


[LON, LAT ] = meshgrid(lon, lat);


% put this 2-month mean SST field into the followint filter to get the
% large scale SST gradient map:
if strcmp(method, 'planefit')
    valid = ~isnan(SST_maps_local_mean);
    C =planefit(LON(valid), LAT(valid), SST_maps_local_mean(valid));
    SST_LS = C(1)*LON + C(2)*LAT + C(3);
    
elseif strcmp(method, 'spectrum')
    % double check if nan exsits:
    if any(isnan(SST_maps_local_mean(:)))
%         figure(1);
%         clf; subplot(1,2,1); imagesc(SST_maps_local_mean);
%     % if so, have to fill it as well...
%        valid = ~isnan(SST_maps_local_mean);
%        F = scatteredInterpolant(LON(valid), LAT(valid), SST_maps_local_mean(valid),'natural','nearest');
%        SST_maps_local_mean(~valid) = F(LON(~valid), LAT(~valid));
%        
%         subplot(1,2,2); imagesc(SST_maps_local_mean); title('after interp.');
       % land mask: 
       % use the closest value:
       valid = ~isnan(SST_maps_local_mean);
       %F = scatteredInterpolant(LON(valid), LAT(valid), SST_maps_local_mean(valid),'natural','nearest');
       SST_maps_local_mean(~valid) = mean(SST_maps_local_mean(valid),'omitnan');
       %SST_maps_local_mean(~valid) = F(LON(~valid), LAT(~valid));
    end
    % compute the 2D spectrum:
    spatial_meanSST = mean(SST_maps_local_mean,'omitnan');
    FT = fft2(SST_maps_local_mean - spatial_meanSST );
    PF = abs(fftshift(FT));
    FT_shifted = fftshift(FT);
    %figure(11); subplot(2,1,1);imagesc(PF);
    
    % compute the wavenumber associated with the map:
    [Ny, Nx] = size(SST_maps_local_mean);
    
    % -- 1.a construct wavenumber space:
    if mod(Nx,2)==0    % even number:
        ix = linspace(-Nx/2,Nx/2-1,Nx);
        %ix = linspace(-Nx/2+1,Nx/2,Nx);
        nxshift = Nx/2;
    else
        ix = linspace(-floor(Nx/2),floor(Nx/2),Nx);
        nxshift = floor(Nx/2)+1;
    end
    
    if mod(Ny,2)==0
        iy = linspace(-Ny/2,Ny/2-1,Ny);
        %iy = linspace(-Ny/2+1,Ny/2,Ny);
        nyshift = Ny/2;
    else
        iy = linspace(-floor(Ny/2),floor(Ny/2),Ny);
        nyshift = floor(Ny/2)+1;
    end
    
    dx = mean(diff(lon))*111*cosd(mean(lat));              % units: km
    dy = mean(diff(lat))*111;
    kx = 1/(Nx*dx) .* ix;        % the last kx is the nyquist frequency, units: cycle per km.
    ky = 1/(Ny*dy) .* iy;
    
    [KX, KY] = meshgrid(kx,ky);
    %
    %figure(11); subplot(2,1,2);imagesc(kx,ky, PF);
    
    % only take variance at scale large to the specified scale.
    k_cutoff = 1/(LS_cutoff);          % cycle per km
    lpf = (KX.^2 + KY.^2)<=k_cutoff.^2;      % the wavenumber magnitude >cutoff will be 0.
    guassian_lpf = 1*exp(-(KX.^2 + KY.^2)./(2*k_cutoff^2));

    % re-construct the SST field using the selected Fourier components: %
    % inverse Fourier transformation.
    %kx_mask = abs(kx)>k_cutoff;
    %ky_mask = abs(ky)>k_cutoff;
    %FT_filtered = FT_shifted;
    %FT_filtered(ky_mask, kx_mask)= 0;        % the way I did it at first give a rectangular mask
    %FT_filtered = FT_shifted.*double(lpf);    % wavenumber magnitude <=cutoff wn will be retained. (circular region)
    FT_filtered = FT_shifted.*guassian_lpf;
    % now, reverse fftshift on the FT_filtered matrix:
    % shift is needed in two dimensions:
    % confirmed with test that the following shift works correctly.
    FT_filtered_shiftback = circshift(FT_filtered, nxshift, 2);      % shift in x dir
    FT_filtered_shiftback = circshift(FT_filtered_shiftback, nyshift, 1);   % shift in y direction
    SST_LS = ifft2(FT_filtered_shiftback) + spatial_meanSST ;

end

SST_LS(isnan(SST_maps_local_mean0))=NaN;

if checkflag
    % make plot to confirm:
    VX = [-60, -48, -48, -60, -60];
    VY = [8, 8, 18, 18,6];
    figure(16);
    subplot(1,3,1)
    pcolor(LON, LAT, SST_maps_local_mean); shading flat; colorbar;
    hold on
    patch(VX',VY','r','faceColor','none');
    axis('square');
    title({'Jan-Feb mean SST map'; '(excluded nan)'});
    caxis([298, 301]);
    
    subplot(1,3,2)
    pcolor(LON, LAT, SST_LS); shading flat; colorbar;
    hold on
    patch(VX',VY','r','faceColor','none');
    
    axis('square');
    
    title({'Jan-Feb mean large scale SST gradient' ;'(from valid data point)'});
    caxis([298, 301]);
    
    subplot(1,3,3)
    pcolor(LON, LAT, SST_maps(:,:,1) - SST_LS); shading flat; colorbar
    hold on
    patch(VX',VY','r','faceColor','none');
    
    axis('square');
    title({' SST anomaly';' at 12UTC'});
    caxis([-1,1]);
end

return