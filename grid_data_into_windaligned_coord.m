function gridded_data = grid_data_into_windaligned_coord(cloud_struct, gridX, gridY)

% inputs: vectors that is already in the wind-aligned coordinate but are
% scatters without a regular coordinate.

nblobs = length(cloud_struct);

method = 'natural';
extrapmethod = 'none';

%gridX = [-3.5:0.1:3.5];
%gridY = [-3.5:0.1:3.5];
[XX, YY] = meshgrid(gridX, gridY);


for i = 1:nblobs
    xA = double(cloud_struct(i).WindAligned_Coord(1,:));
    yA = double(cloud_struct(i).WindAligned_Coord(2,:));
    time = cloud_struct(i).time;
    
    if ~iscolumn(xA)
        xA = xA'; yA= yA'; 
    end
    
    if ~isnan(cloud_struct(i).cloudfreq)
        cloudfreq = cloud_struct(i).cloudfreq;
        cloudycnt = cloudfreq * length(time);
        
        
        if ~iscolumn(cloudfreq)
            cloudfreq = cloudfreq';
        end

        GF = scatteredInterpolant(xA, yA, cloudfreq, method, extrapmethod );
        GF_cc = scatteredInterpolant(xA,yA, cloudycnt, 'nearest','none');
        
        CF_gridded(:,:,i) = GF(XX,YY);
        CloudCnt_gridded(:,:,i) = GF_cc(XX,YY);
        samplesz(i) = length(time);
    else
        
        if ~isempty(cloud_struct(i).SST_cloudmask)
        cloudycnt = double(cloud_struct(i).SST_cloudmask);
        samplesz(i) = length(time);
        CF_gridded(:,:,i) = nan(size(XX));
        else
            cloudycnt = nan(size(XX));
            samplesz(i) = 0;
            CF_gridded(:,:,i) = nan(size(XX));
        end
    end
end

gridded_data.XX = XX;
gridded_data.YY = YY;
gridded_data.cloudyfreq = CF_gridded;
gridded_data.cloudycnt = CloudCnt_gridded;
gridded_data.samplesz = samplesz;


return