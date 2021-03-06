function [composites, sample_size, varargout] = composite_data_in_windaligned_coord_warm_and_cold_updated(dataIn, subset_ids, varargin) %, viewflag, methodstr)
% Purpose: composite the input data (structure) into different subsets
%          using the provided indices:
% Inputs:  dataIn (matlab structure)
%          subset_ids  (required, matlab array)
%          Optional inputs:
%          viewflag : control plotting composite results;
%          methodstr (this is specific for cloudiness computation)
%
% Outputs: composites, sample_size 
%          include the data used to calculate the composite:
%          optional: (coord, hfig), for plotting purposes.
% update: Apr 7, 2022 (no need to regrid data in this function, it has been taken care of before the data exploration.)
%        should I embed the view_composite function here?
%
%        Apr 22: need to update land-sea mask to mask sure that the results
%        arenot affected by land for instance. (I don't think it would
%        since I only select blobs that are identified within the search
%        region..
%
%        June 7, 2022: after updated the format of the input data structure
%        , this function can be simplified a bit. 

%--- 0. input parser --%
%default_method = 'aveCF';     % averaging cloud frequency.
default_flag = 'true';
default_anomaly_type = 'climatological' ;   % 2month, anomaly defined as the difference from the 2-month (or full data period) mean at the feature location.
default_coord.XX = -3:0.05:3;
default_coord.YY = -3:0.05:3;

h = findobj('type','figure');
if ~isempty(h)
    fignum = max([h.Number])+1;
else
    fignum=1;
end
hfig = figure(fignum);
default_ax =gca;


p = inputParser;
addRequired(p, 'dataIn', @isstruct);
addRequired(p, 'subset_ids', @isnumeric);   % can not be a cell array.
addParameter(p, 'viewflag', default_flag, @islogical);
checkString = @(s) any(strcmpi(s, {'aveCF','fromCC'}));
%addParameter(p, 'methodstr', default_method, checkString);
addParameter(p, 'axes', default_ax, @ishandle);
addParameter(p, 'casename', 'none', @ischar);
checkAnomalyStr = @(s) any(strcmpi(s, {'climatological', 'spatial'}));
addParameter(p, 'anomaly_type',default_anomaly_type, checkAnomalyStr);
addParameter(p, 'coord', default_coord, @isstruct);


parse(p, dataIn, subset_ids, varargin{:});
dataIn = p.Results.dataIn;
subset_ids = p.Results.subset_ids;
viewflag = p.Results.viewflag;
%methodstr = p.Results.methodstr;
ax = p.Results.axes;
casename = p.Results.casename;
anomaly_type = p.Results.anomaly_type;
coord = p.Results.coord;

% --------- Start composite (averaging) ------------ %
sample_size.all = length(subset_ids);

% get windaligned grid (it is already contained in the input data
% structure;
% coord.XX = dataIn(1).cloudInfo(1).regridded.XX;
% coord.YY = dataIn(1).cloudInfo(1).regridded.YY;

% use subset_ids to get cloudiness, SST related variables out to do averaging:
%  - aggregate all data together first:
all_cloudInfo = reorganize_struct_to_array([dataIn.cloudInfo]);
all_SSTInfo = reorganize_struct_to_array([dataIn.SSTInfo]);
all_windInfo = reorganize_struct_to_array([dataIn.windInfo]);
all_blobProps = [dataIn.properties];

all_blobmasks_coord = [all_blobProps.blob_image_GeoCoord];
all_blobmasks = cat(2, all_blobProps.blob_image);
all_blobcen = cat(1,all_blobProps.GeoLocs);
tmp = [all_blobProps.stats_selected];
all_blobRad = [tmp.EqvDiam]./2*0.05;

% environmental parameters that would be of interest:
all_maxSSTa = [all_blobProps.max_SSTa];
all_SSTbg = [all_blobProps.ave_SSTbg];
all_blobtime = [];
for i = 1:length(all_blobProps)
    all_blobtime= [all_blobtime, all_blobProps(i).time.*ones(size(all_blobProps(i).max_SSTa))];
end
tmp = [all_blobProps.ave_windInfo];
all_wspd = [tmp.wndspd];
all_wdir = [tmp.wnddir];

% extract cutouts:
all_cloudCOs = [dataIn.cutouts_cloud];
all_windCOs = [dataIn.cutouts_wind];
all_SSTCOs = [dataIn.cutouts_SST];

% the following three needs to be changed.
% establish a nested function to help with subseting the data array;

cloudInfo_sel = get_subset(all_cloudInfo, 3, subset_ids);
SSTInfo_sel = get_subset(all_SSTInfo, 3, subset_ids);
windInfo_sel = get_subset(all_windInfo, 3, subset_ids);
% ---------
blobcoord_sel = all_blobmasks_coord(subset_ids);
blobmask_sel = all_blobmasks(subset_ids);
blobcen_sel = all_blobcen(subset_ids,:);
blobRad_sel = all_blobRad(subset_ids);

cloudCOs_sel = all_cloudCOs(subset_ids);
SSTCOs_sel = all_SSTCOs(subset_ids);
windCOs_sel = all_windCOs(subset_ids);

maxSSTa_sel = all_maxSSTa(subset_ids);
SSTbg_sel = all_SSTbg(subset_ids);
blobtime_sel = all_blobtime(subset_ids);
mwspd_sel = all_wspd(subset_ids);
mwdir_sel = all_wdir(subset_ids);

% indices to select warm and cold features;
%%%% Actually, this is correct.  b/c the gridded data has the same length as subset_ids%%%%%%
Ids.warm = find(all_maxSSTa(subset_ids)>0);
Ids.cold = find(all_maxSSTa(subset_ids)<0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blob_type = fieldnames(Ids); 


% sanity operation: execute only when the blob regridded coordinate is not consistent with
%  the coordinate we wanted it to be on...
% need to input coord_XX...
% if length(coord_XX)~=length(coord.XX )
%     if length(coord_XX)> length(coord.XX)
%         % select subset:
%         xmask = coord_XX(1,:)>=coord.XX(1) & coord_XX(1,:)<=coord.XX(end);
%         ymask = coord_YY(:,1)>=coord.YY(1) & coord_YY(:,1)<=coord.YY(end);
%         cloudInfo_sel = get_subset(cloudInfo_sel, [1,2], {ymask; xmask});
%         windInfo_sel = get_subset(windInfo_sel, [1,2], {ymask; xmask});
%         SSTInfo_sel = get_subset(SSTInfo_sel, [1,2], {ymask; xmask});
%         
%     else
%         return
%     end
% end



% further select the warm and cold blob information out:
for i = 1:2
    BT = blob_type{i};
    sample_size.(BT) = length(Ids.(BT));

    blobs.(BT).cloudInfo = get_subset(cloudInfo_sel, 3, Ids.(BT));
    blobs.(BT).windInfo = get_subset(windInfo_sel, 3, Ids.(BT));
    blobs.(BT).SSTInfo = get_subset(SSTInfo_sel, 3, Ids.(BT));
    
    blobs.(BT).Props.mask = blobmask_sel(Ids.(BT));
    blobs.(BT).Props.mask_coord = blobcoord_sel(Ids.(BT));
    blobs.(BT).Props.blob_cen = blobcen_sel(Ids.(BT),:);
    blobs.(BT).Props.blob_EqvRadDgr  = blobRad_sel(Ids.(BT));
    blobs.(BT).Props.maxSSTa = maxSSTa_sel(Ids.(BT));
    blobs.(BT).Props.SSTbg = SSTbg_sel(Ids.(BT));
    blobs.(BT).Props.time = blobtime_sel(Ids.(BT));
    blobs.(BT).Props.ave_wspd = mwspd_sel(Ids.(BT));
    blobs.(BT).Props.ave_wdir = mwdir_sel(Ids.(BT));
    
    blobs.(BT).cutouts_cldfrac = cloudCOs_sel(Ids.(BT));
    blobs.(BT).cutouts_winddiv = windCOs_sel(Ids.(BT));
    blobs.(BT).cutouts_SSTanom = SSTCOs_sel(Ids.(BT));
end


% construct a 3D matrix from the matlab structure with length of nblobs: (nX, nY, nblobs)
% need to update the field name here:
cloud_vars = {'cloudfrac', 'cloudfrac_highfreq','cloudfrac_relanom_d', ...
              'cloudfrac_relanom_2mon', 'refCF_LP2monmean','refCF_LPdaily'};       % hard coded because I know the data structure. 
SST_vars = {'SST_cutouts', 'SSTa_cutouts', 'eSSTgrad'};
wind_vars = {'winddiv','winddiv_highfreq'};


% the following 3 code blocks are no longer needed.
if 1==0
    for iv = 1:length(cloud_vars)
        CVN = cloud_vars{iv};
        for ib = 1:sample_size.all
            coord_XX = cloudInfo_sel(ib).regridded.XX;
            coord_YY = cloudInfo_sel(ib).regridded.YY;
            
            if length(coord_XX)~=length(coord.XX )
                if length(coord_XX)> length(coord.XX)
                    % select subset:
                    xmask = coord_XX(1,:)>=coord.XX(1) & coord_XX(1,:)<=coord.XX(end);
                    ymask = coord_YY(:,1)>=coord.YY(1) & coord_YY(:,1)<=coord.YY(end);
                    gridded_data.(CVN)(:,:,ib) = cloudInfo_sel(ib).regridded.(CVN)(ymask, xmask);
                    
                else
                    return
                end
                
                
            else
                % this step here can now be removed.
                gridded_data.(CVN)(:,:,ib) = cloudInfo_sel(ib).regridded.(CVN);
            end
        end
    end
    
    for iv = 1:length(wind_vars)
        WVN = wind_vars{iv};
        for ib = 1:sample_size.all
            coord_XX = windInfo_sel(ib).regridded.XX;
            coord_YY = windInfo_sel(ib).regridded.YY;
            
            if length(coord_XX)~=length(coord.XX )
                if length(coord_XX)> length(coord.XX)
                    % select subset:
                    xmask = coord_XX(1,:)>=coord.XX(1) & coord_XX(1,:)<=coord.XX(end);
                    ymask = coord_YY(:,1)>=coord.YY(1) & coord_YY(:,1)<=coord.YY(end);
                    gridded_data.(CVN)(:,:,ib) = windInfo_sel(ib).regridded.(WVN)(ymask, xmask);
                    
                else
                    return
                end
                
                
            else
                gridded_data.(WVN)(:,:,ib) = windInfo_sel(ib).regridded.(WVN);
            end
        end
    end
    
    for iv = 1:length(SST_vars)
        SVN = SST_vars{iv};
        for ib = 1:sample_size.all
            coord_XX = SSTInfo_sel(ib).regridded.XX;
            coord_YY = SSTInfo_sel(ib).regridded.YY;
            
            if ~isempty(SSTInfo_sel(ib).regridded)
                
                if length(coord_XX)~=length(coord.XX )
                    if length(coord_XX)> length(coord.XX)
                        % select subset:
                        xmask = coord_XX(1,:)>=coord.XX(1) & coord_XX(1,:)<=coord.XX(end);
                        ymask = coord_YY(:,1)>=coord.YY(1) & coord_YY(:,1)<=coord.YY(end);
                        gridded_data.(CVN)(:,:,ib) = SSTInfo_sel(ib).regridded.(SVN)(ymask, xmask);
                        
                    else
                        return
                    end
                    
                    
                else
                    gridded_data.(SVN)(:,:,ib) = SSTInfo_sel(ib).regridded.(SVN);
                end
                
            else
                % nan the entry;
                gridded_data.(SVN)(:,:,ib) = gridded_data.(SVN)(:,:,ib-1);
                gridded_data.(SVN)(:,:,ib) = NaN;
                
            end
        end
    end
end

%% ready to compute composite for warm and cold features separately.
% - for clouds:
% I can do the composite directly from the structure I generated above. 
% mean(blobs.(BT).cloudInfo.(FN),'');
% median(blobs.(NT).cloudInfo.(FN), '');
% clean the original code up a bit better.

% compute the mean for these different structure. 
timedimID = 3;
for i = 1:2
    BT = blob_type{i};
    cldfrac_comp.(BT) = get_statistics_from_struct(blobs.(BT).cloudInfo, 'mean', timedimID);
    winddiv_comp.(BT) = get_statistics_from_struct(blobs.(BT).windInfo, 'mean', timedimID);
    SST_comp.(BT) = get_statistics_from_struct(blobs.(BT).SSTInfo, 'mean', timedimID);    
end

% combine the fields from these three data together:
cldvars = fieldnames(cldfrac_comp.warm);
sstvars = fieldnames(SST_comp.warm);
windvars = fieldnames(winddiv_comp.warm);
for i = 1:2
    BT = blob_type{i};
    for iv = 1:length(cldvars)
        VN = cldvars{iv};
        composites.(BT).(VN) = cldfrac_comp.(BT).(VN);
    end
    
    for iv = 1:length(sstvars)
        VN = sstvars{iv};
        composites.(BT).(VN) = SST_comp.(BT).(VN);
    end
    
    for iv = 1:length(windvars)
        VN = windvars{iv};
        composites.(BT).(VN) = winddiv_comp.(BT).(VN);
    end
end




%% now, do some plotting in the following format:
% SSTa , effecitve SST grad, cloudiness, cloudiness enhancement:
% for cold and warm SST; 
if viewflag
    % 
    cmap_div = redblue(32);
    cmap_jet = jet(32);
    varlist = {'SSTa_cutouts','eSSTgrad', 'cloudfrac', ...
                'cloudfrac_relanom_2mon','winddiv_highfreq'}; %,'cloudfreq_std', 'cldfreq_anom_std'};     % why I still don't have cloudfreq_anom, needs to be added.??
    varname = {'SST anom.', '\bf{u} \cdot \nablaSST', 'cloudiness', ...
               {'cloudiness'; 'spatial anomalies'}, {'wind divergence'}}; %{'cloudiness anomalies', 'from "climatology"' }%, 'cloudiness stdv.', 'cloudiness anomaly stdv.'};
    varunits = {'K', 'K/hr', '', '', 's^{-1}'};
    colormaps = {cmap_div, cmap_div, cmap_jet, cmap_div, cmap_div};
    caxis_range_default ={[-0.15, 0.15], [-7.5*10^(-5), 7.5*10^-5], [0.3, 0.5], [-0.1,0.1], [-5e-4, 5e-4]};    %[-0.3,0.3]
    
    
    nrow = 2;  ncol=length(varlist);
    xspace = 0.08; yspace = 0.08;
    pos = customize_subplot_size(nrow, ncol, xspace, yspace);
    
    figure(hfig.Number);
    set(hfig, 'pos',[339, 410, 1280,600]);
    set(hfig, 'Name', casename);
    for ir = 1:nrow
        BT = blob_type{ir};
        
        for icol = 1:ncol      % warm and cold features
            VN = varlist{icol};
            val = composites.(BT).(VN);
            maxval = max(abs(val(:)));
            caxis_range = [-maxval, maxval];
            
            ip = icol+ncol*(ir-1);
            hsub = subplot(nrow, ncol, ip);   % different variables:

            pcolor(coord.XX, coord.YY, val); shading flat;
            hold on
%             if icol == 3
%                 if any(strcmpi(casename, {'flowers', 'flower'}))
%                     cldlev=[0.5:0.02:0.7];
%                 else
%                     cldlev = [0.3:0.02:0.5];
%                 end
%                 
%                 %[c,hl] = contour(coord.XX, coord.YY, val, cldlev,'k');
%                 %clabel(c,hl, cldlev,'fontsize',6,'LabelSpacing',600);
%             end
%             if icol >3
%                 if any(strcmpi(casename, {'flowers', 'flower'}))
%                     clevs = [-0.3:0.05:0.3]; %[-1.5:0.1:1.5];
%                 else
%                     clevs = [-0.3:0.05:0.3];
%                 end
%                 %[c,hl] = contour(coord.XX, coord.YY, val, clevs,'k');
%                 %clabel(c,hl, clevs,'fontsize',6,'LabelSpacing',600);
%             end
            hold on
            circle(0,0,1,'k',2);
%             if contains(VN, {'cloudfrac'}); %,'cldfreq'}) 
%                 cmap = cmap_jet;
%             else
%                 if icol<=4
%                     cmap = cmap_div;
%                 else
%                     cmap = cmap_jet;
%                 end
%             end
            cmap = colormaps{icol};
            colormap(hsub, cmap);
            hb = colorbar(hsub); 
            set(get(hb,'xlabel'),'String', varunits{icol});
            caxis(caxis_range_default{icol});
            
%             if icol>3 & icol<5
%                 if any(strcmpi(casename, {'flowers', 'flower'}))
%                     set(hb,'xtick',[-0.5:0.1:0.5]); %[-1.5:0.3:1.5]);
%                 else
%                     set(hb,'xtick',[-0.5:0.1:0.5]);
%                 end
%             end
%             
%             if icol ==5
%                 set(hb,'xtick',[-0.5:0.1:0.5].*10^(-2));
%             end
%             


            if mod(ip,ncol)==1 
                titlestr = {varname{icol} ; ...
                          ['samples (' num2str(sample_size.(BT)) ')']};
            else
                titlestr = varname{icol};
            end
            
            title(titlestr,'fontsize',10.5);
%             if icol<3
%                 caxis(caxis_range);
%             else
%                 if icol==3
%                     if any(strcmpi(casename, {'flowers', 'flower'}))
%                         caxis([0.5, 0.7]);
%                     else
%                         caxis([0.3 0.5]);
%                     end
%                 elseif icol>3 & icol<5
%                     %if icol<=4
%                         if any(strcmpi(casename, {'flowers', 'flower'}))
%                             caxis([-0.3, 0.3]);
%                             %caxis([-1.5,1.5]);
%                         else
%                             caxis([-0.3, 0.3]);
%                         end
%                 else %icol==5
%                     caxis(caxis_range_default{5});
%                     
%                 end
%             end
            axis('square');
            if mod(ip, ncol)==1
                ylabel({'downwind'; 'normalized distance'});
            end
            
            if ceil(ip/ncol)==nrow
                xlabel({'crosswind'; 'normalized distance'});
            end
            set(gca,'fontsize',10);
            
            set(hsub,'pos', pos{ip});
        end
    end
    
else
    close(hfig)
end


% output:
varargout{1} = coord;
varargout{2} = blobs;
if viewflag
    varargout{3} = hfig;
end
       

end

% --------- util function: ------------ %
function struct_subseted = get_subset(structIn, dimN, subset_ids)

fieldn = fieldnames(structIn);
for iv = 1:length(fieldn)
    FN = fieldn{iv};
    
    if length(dimN)==1
        if dimN==3
            struct_subseted.(FN) = structIn.(FN)(:,:,subset_ids);
        elseif dimN == 2
            struct_subseted.(FN) = structIn.(FN)(:,subset_ids,:);
        elseif dimN ==1
            struct_subseted.(FN) = structIn.(FN)(subset_ids,:,:);
        end
    
    elseif length(dimN)==2
        if dimN(1)==1 && dimN(2) ==2
            dim1_ids = subset_ids{1};
            dim2_ids = subset_ids{2};
            struct_subseted.(FN) = structIn.(FN)(dim1_ids,dim2_ids,:);
        end
    end
        
    
    
end

end

function structOut = get_statistics_from_struct(structIn, func, dimN)
fieldn = fieldnames(structIn);
for iv = 1:length(fieldn)
    FN = fieldn{iv};
    
    structOut.(FN)= feval(func, structIn.(FN), dimN, 'omitnan');
    
end

end

function mean_tmp = mean_omitnan_with_threshold(tmp_data, dimID, nanprc)
nt = size(tmp_data,3);
nanmask = isnan(tmp_data);
% check how much data is missing at each grid point:
missing_data_nt = sum(double(nanmask),dimID)./nt;
thres = missing_data_nt< nanprc;
thres_3D = repmat(thres,1,1,nt);
figure(10);clf;
pcolor(double(thres)); colorbar;
pcolor(missing_data_nt); shading flat; colorbar;

mean_tmp = nan(size(thres));

masked_datavec = tmp_data(thres_3D);
masked_data3D = reshape(masked_datavec, size(thres,1), size(thres,2), nt);

mean_tmp(thres)  = mean(masked_data3D, dimID,'omitnan');

end


% if 1==0
%     % Obsolete:
%     for ic = 1:2
%         BT = blob_type{ic};
%         m =Ids.(BT);
%         sample_size.(BT) = length(m);
%         
%         if strcmp(methodstr,'aveCF')
%             
%             %CF_tmp = mean_omitnan_with_threshold(gridded_data.cloudfrac(:,:,m), 3, nanprc);
%             CF_tmp  = mean(gridded_data.cloudfrac(:,:,m), 3,'omitnan');
%             
%         elseif strcmp(methodstr, 'fromCC')
%             CF_tmp = sum(gridded_data.cloudcnt(:,:,m), 3)./sample_size.(BT);
%             
%         end
%         CF_std = std(gridded_data.cloudfrac(:,:,m),1,3,'omitnan');
%         
%         % mask out land: (always ==1); (not correct here, because the blobs are
%         % not at the same location, this criterion applies in the 2D datasets.
%         % )
%         % CF_tmp(CF_tmp==1)=NaN;
%         % CF_var(CF_var==0) =NaN;
%         composites.cloudfrac.(BT) = CF_tmp;
%         composites.cloudfrac_std.(BT) = CF_std;
%         composites.cloudfrac_samples.(BT) = gridded_data.cloudfrac(:,:,m);
%         
%         %if strcmpi(anomaly_type,'climatological')
%         
%         CFa_clim = mean(gridded_data.cloudfrac_relanom_2mon(:,:,m), 3, 'omitnan');
%         %CFa_clim = mean_omitnan_with_threshold(gridded_data.cldfrac_anom(:,:,m), 3, nanprc);
%         %% I can also NaN out results where NaN exist more than certain percentage of the data;
%         
%         %CFa_tmp(isnan(CF_tmp))=NaN;    % mask out
%         CFa_clim_std = std(gridded_data.cloudfrac_relanom_2mon(:,:,m), 1,3,'omitnan');
%         
%         % elseif strcmpi(anomaly_type, 'spatial')
%         spatial_CF_means =  mean(gridded_data.cloudfrac(:,:,m), [1,2],'omitnan');
%         spatial_CF_means_3D = repmat(spatial_CF_means, size(gridded_data.cloudfrac,1), size(gridded_data.cloudfrac,2),1);
%         CFa_tmp_all = (gridded_data.cloudfrac(:,:,m) -spatial_CF_means_3D) ./spatial_CF_means_3D;
%         CFa_spatial_tmp = mean(CFa_tmp_all, 3,'omitnan');
%         %CFa_spatial_tmp = mean_omitnan_with_threshold(CFa_tmp_all, 3, nanprc);
%         %CFa_tmp(isnan(CF_tmp))=NaN;
%         CFa_spatial_std = std(CFa_tmp_all, 1, 3,'omitnan');
%         %end
%         
%         composites.cldfrac_anomSpa.(BT) = CFa_spatial_tmp;
%         composites.cldfrac_anomSpa_std.(BT) = CFa_spatial_std;
%         
%         composites.cldfrac_anomClim.(BT) = CFa_clim;
%         composites.cldfrac_anomClim_std.(BT) = CFa_clim_std;
%         
%         composites.cldfracAnom_samples.(BT) = gridded_data.cloudfrac_relanom_2mon(:,:,m);
%         
%         
%         % if exist refCF_xxx field:
%         if any(contains(fieldnames(gridded_data), 'refCF'))
%             for iv = 3:4
%                 VN = cloud_vars{iv};
%                 composites.(VN).(BT) = gridded_data.(VN)(:,:,m);
%             end
%         end
%         
%         % - for SSTs:
%         for iv = 1:length(SST_vars)
%             SVN = SST_vars{iv};
%             composites.(SVN).(BT) = mean(gridded_data.(SVN)(:,:,m), 3,'omitnan');  %  omitted NaN, since those will be land
%             %composites.(SVN).(BT) = mean_omitnan_with_threshold(gridded_data.(SVN)(:,:,m), 3, nanprc);
%             SVN_samples = [SVN '_samples'];
%             composites.(SVN_samples).(BT) = gridded_data.(SVN)(:,:,m);
%         end
%         
%         % - for wind divergence:
%         for iv = 1:length(wind_vars)
%             WVN = wind_vars{iv};
%             composites.(WVN).(BT) = mean(gridded_data.(WVN)(:,:,m),3,'omitnan');
%             WVN_samples = [WVN '_samples'];
%             composites.(WVN_samples).(BT) = gridded_data.(WVN)(:,:,m);
%         end
%         
%         
%         %     %% This may not be the right place to do this calculation.
%         %     %% significant level estimate: it is a function of x and y, and this can be done for each parameter (esp. cluodiness and wind div.)
%         %     % use bootstrapping again to do this;
%         if 1==0
%             samplesz_prc = 10;     % 10 percent of data;
%             sigalpha = 0.05;
%             nrepeat = 1000;
%             SigLev2D = estimate_significant_level_from_subsets(gridded_data.(WVN)(:,:,m), samplesz_prc, sigalpha, nrepeat);
%             SigLev2D_cloud = estimate_significant_level_from_subsets(gridded_data.cloudfrac_relanom_2mon(:,:,m), samplesz_prc, sigalpha, nrepeat);
%         end
%         
%     end
% end