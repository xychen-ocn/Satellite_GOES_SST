function gridded_data = grid_data_into_windaligned_coord(data_struct, gridX, gridY)

% inputs: vectors that is already in the wind-aligned coordinate but are
% scatters without a regular coordinate.
% Dev Notes:
%  - Apr 6, 2022: generalize function to deal with different data (e.g.,
%  cloudiness, SST, wind), will need to understand the data fields in different dataset a bit better 
% look at Park et al. (2006) again to see how they deal with change of coordinate for cloudiness.
%  - In Park, the values are averaged within each normalized distance bin.
%  (will need to write a different function to do this procedure myself.)
%  - I should probably write something similar to double check if that
%  changes the results at all.

nblobs = length(data_struct);

method = 'natural';
extrapmethod = 'none';

%gridX = [-3.5:0.1:3.5];
%gridY = [-3.5:0.1:3.5];
[XX, YY] = meshgrid(gridX, gridY);

% need to make sure that cloud related variables are reshaped into vector
% before doing interpolations. 


for i = 1:nblobs
    xA = double(data_struct(i).WindAligned_Coord(1,:));
    yA = double(data_struct(i).WindAligned_Coord(2,:));
    
    if ~iscolumn(xA)
        xA = xA'; yA= yA';
    end
    
    if isfield(data_struct,'cloudfreq') 
        time = data_struct(i).time;
        
        % deal with cloudfreq; (obsolete..)
        if ~isnan(data_struct(i).cloudfreq)
            cloudfreq = reshape(data_struct(i).cloudfreq,[],1);
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
        
        
        % deal with cloudfreq_anom:  %
        %% red flag up: I might have forgotten to add (i) to the following line.. (shouldn't caused problem because data_struct is storing information for only 1 blob when this functino is called...
        dataval = reshape(data_struct(i).cloudfreq_anom, [],1);   % should I add (i) here? (double check is required, this may results to the "wrong" results..)
        if ~isnan(dataval)
            F_varn = scatteredInterpolant(xA, yA, dataval, method, extrapmethod);
            gridded_data.cldfreq_anom(:,:,i) = F_varn(XX,YY);
        else
            gridded_data.cldfreq_anom(:,:,i) = nan(size(XX));
        end
        
        
        
    elseif isfield(data_struct, 'SST_cutouts')  || isfield(data_struct, 'winddiv')
        % do operation for two variables in this data structure;
        if isfield(data_struct, 'SST_cutouts')
            var_list = {'SST_cutouts','SSTa_cutouts', 'eSSTgrad','refCF_LP2monmean','refCF_LPdaily'};
        else
            var_list = {'winddiv','winddiv_highfreq', 'uwnd','vwnd'};
        end
        for iv = 1:length(var_list)
            % directly do interpolation in a similar way without too much
            % concern as above.
            if isfield(data_struct, var_list{iv})
                dataval = data_struct(i).(var_list{iv});
                if numel(dataval) ~= length(dataval)
                    dataval = reshape(dataval, size(xA));
                else
                    if ~iscolumn(dataval)
                        dataval = dataval';
                    end
                end
            
                
                if length(dataval)~=1 & ~isnan(dataval)
                    F_varn = scatteredInterpolant(xA, yA, dataval, 'nearest', extrapmethod);
                    gridded_data.(var_list{iv})(:,:,i) = F_varn(XX,YY);
                else
                    gridded_data.(var_list{iv})(:,:,i) = nan(size(XX));
                end
            end
            
        end
        
    elseif isfield(data_struct, 'cloudfrac')
        
        % provide a variable list here to regridd all the variables in the
        % data structure. 
        var_list = {'cloudfrac','cloudfrac_highfreq', 'cloudfrac_relanom_d', ...
                    'cloudfrac_relanom_2mon', 'refCF_LP2monmean', 'refCF_LPdaily'};
        % do the interpolation on to the normalized grid using nearest
        % method again:
         for iv = 1:length(var_list)
                % directly do interpolation in a similar way without too much
                % concern as above.
                if isfield(data_struct, var_list{iv})
                    dataval = double(data_struct(i).(var_list{iv}));
                    if ~isempty(dataval) & all(~isnan(dataval))
                        
                        if numel(dataval) ~= length(dataval)
                            dataval = reshape(dataval, size(xA));
                        else
                            if ~iscolumn(dataval)
                                dataval = dataval';
                            end
                        end
                        
                        
                        
                        F_varn = scatteredInterpolant(xA, yA, dataval, 'nearest', extrapmethod);
                        gridded_data.(var_list{iv})(:,:,i) = F_varn(XX,YY);
                    else
                        gridded_data.(var_list{iv})(:,:,i) = nan(size(XX));
                    end
                end
                
         end
           
        
        
%         cldfrac = data_struct(i).cloudfrac;
%         if all(~isnan(cldfrac)) & ~isempty(cldfrac)
%             if numel(cldfrac) ~= length(cldfrac)
%                 cldfrac = reshape(cldfrac, [],1);
%             else
%                 if ~iscolumn(cldfrac)
%                     cldfrac = cldfrac';
%                 end
%             end
%             F = scatteredInterpolant(xA,yA, cldfrac, 'nearest', extrapmethod);   % extrapolation to none
%             gridded_data.cloudfrac(:,:,i) = F(XX,YY);
%         else
%             gridded_data.cloudfrac(:,:,i) = nan(size(XX));
%         end
%         gridded_data.samplesz = data_struct(i).num_of_hrly_maps;
% 
%         
%         cldfrac_anom = data_struct(i).cloudfrac_anom;
%         if all(~isnan(cldfrac_anom)) & ~isempty(cldfrac)
%             if numel(cldfrac_anom) ~= length(cldfrac_anom)
%                 cldfrac_anom = reshape(cldfrac_anom, [],1);
%             else
%                 if ~iscolumn(cldfrac_anom)
%                     cldfrac_anom = cldfrac_anom';
%                 end
%             end
%             
%             Fa = scatteredInterpolant(xA,yA, cldfrac_anom, 'nearest', extrapmethod);
%             gridded_data.cldfrac_anom(:,:,i) = Fa(XX,YY);
%         else
%             gridded_data.cldfrac_anom(:,:,i) = nan(size(XX));
%         end
%         
%         if any(contains(fieldnames(data_struct), 'refCF'))
%             var_list = {'refCF_LP2monmean','refCF_LPdaily'};
%             for iv = 1:length(var_list)
%                 % directly do interpolation in a similar way without too much
%                 % concern as above.
%                 if isfield(data_struct, var_list{iv})
%                     dataval = data_struct(i).(var_list{iv});
%                     if numel(dataval) ~= length(dataval)
%                         dataval = reshape(dataval, size(xA));
%                     else
%                         if ~iscolumn(dataval)
%                             dataval = dataval';
%                         end
%                     end
%                     
%                     
%                     if length(dataval)~=1 & ~isnan(dataval)
%                         F_varn = scatteredInterpolant(xA, yA, dataval, 'nearest', extrapmethod);
%                         gridded_data.(var_list{iv})(:,:,i) = F_varn(XX,YY);
%                     else
%                         gridded_data.(var_list{iv})(:,:,i) = nan(size(XX));
%                     end
%                 end
%                 
%             end
%             
%         end
            
        
        
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
