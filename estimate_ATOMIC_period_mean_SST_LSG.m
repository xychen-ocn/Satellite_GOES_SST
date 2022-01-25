function [SST_LS, SST_maps_local_mean] = estimate_ATOMIC_period_mean_SST_LSG(varargin)
% purpose: utilize 2 month (Jan, Feb) GOES-16 SST to estimate the large
% scale SST gradient in the field of interest;
%  - 2 methods will be implemented in this function: 
%    1) a simple plane fit;
%    2) spectral reduction, needs to specify a scale to define large scale.

default_method = 'planefit';
default_checkflag = true;
default_LScutoff = 1000;    % km
%% build up the parse object to read in element;

p = inputParser;
% define parsing rules:
checkString = @(s) any(strcmpi(s, {'planefit','spectrum'}));
addParameter(p,'method', default_method, checkString);
addParameter(p, 'CutoffScale', default_LScutoff, @isnumeric);           % units: km
addParameter(p, 'checkflag', default_checkflag, @islogical);

% ready to parse the inputs:
parse(p, varargin{:});
method = p.Results.method;
LS_cutoff = p.Results.CutoffScale;
checkflag = p.Results.checkflag;

%% read all the data:
svdataFN = 'ATOMIC_JanFeb_SSTmaps.mat';
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';

if ~exist([datadir filesep svdataFN],'file') 
    files = dir([datadir filesep 'GOES_SST*.mat']);
    filenames = {files.name};
    for j = 1:length(filenames)
        matFN = filenames{j};
        tmp = strsplit(matFN,'_');
        dateID = tmp{end}(1:5);
        
        load([datadir filesep matFN]);
        
        disp(['--> working on ' dateID]);
        
        % take subset to exclude islands...
        xstart = -59;
        % longitude mask below actually applies to all the time instance:
        lonmask = GOES_ATOMIC.lon(:,1)>xstart;
        datasub.lon = double(GOES_ATOMIC.lon(lonmask,1));
        datasub.lat = double(GOES_ATOMIC.lat(:,1));
        
        [LON, LAT]= meshgrid(datasub.lon, datasub.lat);
        
        % stack up data:
        SST_maps = permute(GOES_ATOMIC.sea_surface_temperature(lonmask,:,:), [2,1,3]);
        nt = size(SST_maps,3);
        
        if j==1
            it0 = 1;
        end
        itN = it0+nt-1;
        SST_all(:,:,it0:itN) = SST_maps;
        time_all(it0:itN) = GOES_ATOMIC.time_num;
        it0 = itN + 1;
    end
    save([datadir filesep svdataFN],'SST_all','time_all', 'LON','LAT','datasub','-v7.3');
    
else
    load([datadir filesep svdataFN]);
end
%% compute mean from valide data at each grid point for the Jan-Feb period.
% the order above is not ascending. But it is ok for the mean field
% calculation.
SST_maps_local_mean = mean(SST_all, 3, 'omitnan');

% put this 2-month mean SST field into the followint filter to get the
% large scale SST gradient map:
if strcmp(method, 'planefit')
    valid = ~isnan(SST_maps_local_mean);
    C =planefit(LON(valid), LAT(valid), SST_maps_local_mean(valid));
    SST_LS = C(1)*LON + C(2)*LAT + C(3);
    
elseif strcmp(method, 'spectrum')
    % compute the 2D spectrum:
    FT = fft2(SST_maps_local_mean);
    %PF = abs(fftshift(FT));
    FT_shifted = fftshift(FT);
    %figure(11); subplot(2,1,1);imagesc(PF);
    
    % compute the wavenumber associated with the map:
    [Ny, Nx] = size(SST_maps_local_mean);
    
    % -- 1.a construct wavenumber space:
    if mod(Nx,2)==0    % even number:
        ix = linspace(-Nx/2,Nx/2-1,Nx);
        nxshift = Nx/2;
    else
        ix = linspace(-floor(Nx/2),floor(Nx/2),Nx);
        nxshift = floor(Nx/2)+1;
    end
    
    if mod(Ny,2)==0
        iy = linspace(-Ny/2,Ny/2-1,Ny);
        nyshift = Ny/2;
    else
        iy = linspace(-floor(Ny/2),floor(Ny/2),Ny);
        nyshift = floor(Ny/2)+1;
    end
    
    dx = mean(diff(datasub.lon))*111*cosd(mean(datasub.lat));              % units: km
    dy = mean(diff(datasub.lat))*111;
    kx = 1/(Nx*dx) .* ix;        % the last kx is the nyquist frequency, units: cycle per km.
    ky = 1/(Ny*dy) .* iy;
    
    %[KX, KY] = meshgrid(kx,ky);
    %figure(11); subplot(2,1,2);imagesc(kx,ky, PF);
    
    % only take variance at scale large to the specified scale.
    k_cutoff = 1/(LS_cutoff);          % cycle per km
    
    % re-construct the SST field using the selected Fourier components: %
    % inverse Fourier transformation.
    kx_mask = abs(kx)>k_cutoff;
    ky_mask = abs(ky)>k_cutoff;
    FT_filtered = FT_shifted;
    FT_filtered(ky_mask, kx_mask)= 0;
    
    % now, reverse fftshift on the FT_filtered matrix:
    % shift is needed in two dimensions:
    % confirmed with test that the following shift works correctly.
    FT_filtered_shiftback = circshift(FT_filtered, nxshift, 2);      % shift in x dir
    FT_filtered_shiftback = circshift(FT_filtered_shiftback, nyshift, 1);   % shift in y direction
    SST_LS = ifft2(FT_filtered_shiftback);

end


% make plot to confirm:

figure();
subplot(1,3,1)
pcolor(LON, LAT, SST_maps_local_mean); shading flat; colorbar;
axis('square');
title({'Jan-Feb mean SST map'; '(excluded nan)'});
caxis([298, 301]);

subplot(1,3,2)
pcolor(LON, LAT, SST_LS); shading flat; colorbar;
axis('square');
title({'Jan-Feb mean large scale SST gradient' ;'(from valid data point)'});
caxis([298, 301]);

subplot(1,3,3)
pcolor(LON, LAT, SST_all(:,:,(28+8)*24+12) - SST_LS); shading flat; colorbar
axis('square');
title({' SST anomaly';' at 12UTC'});
caxis([-1,1]);

return