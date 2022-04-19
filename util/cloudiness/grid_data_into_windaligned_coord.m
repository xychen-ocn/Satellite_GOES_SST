function gridded_data = grid_data_into_windaligned_coord(data_struct, gridX, gridY)

% inputs: vectors that is already in the wind-aligned coordinate but are
% scatters without a regular coordinate.
% Dev Notes:
%  - Apr 6, 2022: generalize function to deal with different data (e.g.,
%  cloudiness, SST, wind), will need to understand the data fields in different dataset a bit better 


nblobs = length(data_struct);

method = 'natural';
extrapmethod = 'none';

%gridX = [-3.5:0.1:3.5];
%gridY = [-3.5:0.1:3.5];
[XX, YY] = meshgrid(gridX, gridY);


for i = 1:nblobs
    xA = double(data_struct(i).WindAligned_Coord(1,:));
    yA = double(data_struct(i).WindAligned_Coord(2,:));
    
    if ~iscolumn(xA)
        xA = xA'; yA= yA';
    end
    
    if isfield(data_struct,'cloudfreq')
        time = data_struct(i).time;
        
        % deal with cloudfreq;
        if ~isnan(data_struct(i).cloudfreq)
            cloudfreq = data_struct(i).cloudfreq;
            cloudycnt = cloudfreq * length(time);
            
            if ~iscolumn(cloudfreq)
                cloudfreq = cloudfreq';
            end
            
            GF = scatteredInterpolant(xA, yA, cloudfreq, method, extrapmethod );
            CF_gridded(:,:,i) = GF(XX,YY);
            
            % is this a good approach?? I think this only make sense if the native grid is coarser than the new grid and not the other way around.
            GF_cc = scatteredInterpolant(xA,yA, cloudycnt, 'nearest','none');
            CloudCnt_gridded(:,:,i) = GF_cc(XX,YY);
            samplesz(i) = length(time);
            
        else   % cloudfreq input is nan due to not enough cloud mask in this day.
            
            if ~isempty(data_struct(i).SST_cloudmask)    % have cloud mask;
                % SST cloud mask is not a logical variable (1/0)
                logical_mask = isnan(data_struct(i).SST_cloudmask);
                cloudycnt = sum(double(logical_mask),2);
                GF_cc = scatteredInterpolant(xA,yA, cloudycnt, 'nearest','none');
                CloudCnt_gridded = GF_cc(XX,YY);
                
                samplesz(i) = length(time);
                CF_gridded(:,:,i) = nan(size(XX));
                
            else    % not cloud mask.
                CloudCnt_gridded(:,:,i) = nan(size(XX));
                samplesz(i) = 0;
                CF_gridded(:,:,i) = nan(size(XX));
            end
        end
        
        
        % deal with cloudfreq_anom:
        dataval = data_struct.cloudfreq_anom;
        if ~isnan(dataval)
            F_varn = scatteredInterpolant(xA, yA, dataval, method, extrapmethod);
            gridded_data.cldfreq_anom(:,:,i) = F_varn(XX,YY);
        else
            gridded_data.cldfreq_anom(:,:,i) = nan(size(XX));
        end
        
        
        
    elseif isfield(data_struct, 'SST_cutouts') 
        % do operation for two variables in this data structure;
        var_list = {'SST_cutouts','SSTa_cutouts', 'eSSTgrad'};
        for iv = 1:length(var_list)
            % directly do interpolation in a similar way without too much
            % concern as above.
            if isfield(data_struct, var_list{iv})
                dataval = data_struct.(var_list{iv});
                if numel(dataval) ~= length(dataval)
                    dataval = reshape(dataval, size(xA));
                end
                F_varn = scatteredInterpolant(xA, yA, dataval, 'nearest', extrapmethod);
                gridded_data.(var_list{iv})(:,:,i) = F_varn(XX,YY);
            end
            
        end
        
    else
        disp('need to be build first, sorry');
        return
    end
end      % end looping through all features.

gridded_data.XX = XX;
gridded_data.YY = YY;
if isfield(data_struct, 'cloudfreq')
    gridded_data.cloudfreq = CF_gridded;
    gridded_data.cloudcnt = CloudCnt_gridded;
    gridded_data.samplesz = samplesz;
end

return