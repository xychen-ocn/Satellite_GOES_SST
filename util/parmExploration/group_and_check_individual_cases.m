function [selected_frac, QInfo] = group_and_check_individual_cases(coord, data,  figsvdir)

% blob_type = 'warm';
% data = blobInfo.(blob_type);
cloudInfo = data.cloudInfo;  % 1 = warm; 2 = cold;
windInfo = data.windInfo;
SSTInfo = data.SSTInfo;
cutouts_CFa = data.cutouts_cldfrac;
cutouts_winddiv = data.cutouts_winddiv;
cutouts_SSTa = data.cutouts_SSTanom;
blobprops = data.Props;
blobprops.blob_cen = blobprops.blob_cen';
blobID = 1:length(cutouts_SSTa);

varlist = {'cloudfrac_highfreq','winddiv_highfreq'};
coordR = sqrt(coord.XX.^2 + coord.YY.^2);

nt = size(cloudInfo.cloudfrac,3);
%coordR3D = repmat(coordR,1,1,nt);
for iv = 1:2
    VN = varlist{iv};
    cenmask = coordR<=0.8;
    for it = 1:nt
        if contains(VN, 'cloud')
            val = cloudInfo.(VN)(:,:,it);
        else
            val = -windInfo.(VN)(:,:,it);
        end
        cencum.(VN)(it) = sum(val(cenmask),'all','omitnan');
    end
end

%% 1. group or cluster different cases;
%   - 2 or 4 groups: within the center of 1R
%     a1) + cldfrac & - windidv; a2) -cldfrac & +windiv;
%     b1) + cldfrac & + winddiv; b2) -cldfrac & - winddiv  

CFa_norm = (cencum.cloudfrac_highfreq-0)./std(cencum.cloudfrac_highfreq);
Conv_norm = (cencum.winddiv_highfreq-0)./std(cencum.winddiv_highfreq);

rthres = 1;  % select cases that have both parameters vary more than 1 standard deviation.
pbound = 0.5;   % parameter selection boundary when r> rthres;

% x_norm = [-round(max(CFa_norm)):0.05:round(max(CFa_norm))];
% y_norm = [-round(max(Conv_norm)):0.05:round(max(Conv_norm))];
r_norm = sqrt(CFa_norm.^2+Conv_norm.^2);
selmask = r_norm >rthres & (abs(CFa_norm)>pbound & abs(Conv_norm)>pbound);


figure(10); clf;
plot(CFa_norm, Conv_norm,'.b');
grid on
xlabel('normalized cloudfrac anormaly');
ylabel('normalized wind convergence anormaly');
hold on
plot([0 0], [-4 4],'--k');
plot([-4 4],[0 0], '--k');
hold on
circle(0,0, rthres, 'r',2);
plot(CFa_norm(selmask), Conv_norm(selmask),'or');

% sort data by its distance to the center, select data points that has
% variation larger than 1 std on both parameter.
% (0,0): means that cumulative sum of cloudfraction anormaly and
% convergence are both 0 with 1R 
% (+,+): means that increase cloud fraction and wind convergence on average
%        within 1R  (-,-) opposite pairs.
% (+,-): means that increase cloud fraction and wind divergence on average
% within 1R; (-,+) opposite pairs.
% (->0, +/-): compensating change in the cloud frac anom, but not in the
% wind pattern.

% now for each quardrant, I can sort this metric by its absolute value
% (distance from 0) and plot the images out;
selected_frac = length(find(selmask==1))/nt * 100;

%% under construction.
% use selmask to get the blobs of interest out:
cloudInfo_sub = get_subset(cloudInfo, 3, selmask);
windInfo_sub = get_subset(windInfo, 3, selmask);
SSTInfo_sub = get_subset(SSTInfo, 3, selmask);
blobprops_sub = get_subset(blobprops, 2, selmask);

cutouts_sub.CFa = cutouts_CFa(selmask);
cutouts_sub.winddiv = cutouts_winddiv(selmask);
cutouts_sub.SSTa = cutouts_SSTa(selmask);
COvarNames = fieldnames(cutouts_sub);

r_norm_sel = r_norm(selmask);
CFa_norm_sel = CFa_norm(selmask);
Conv_norm_sel = Conv_norm(selmask);
blobID_sel = blobID(selmask);

% further select subset to divide data into 4 groups:
quadrant_angle = atan2(Conv_norm(selmask), CFa_norm(selmask))*180/pi;

Qmask{1} = quadrant_angle>0 & quadrant_angle<90;
Qmask{2} = quadrant_angle>90 & quadrant_angle<180;
Qmask{3} = quadrant_angle>-180 & quadrant_angle<-90;
Qmask{4} = quadrant_angle<0 & quadrant_angle>-90;

NQ = length(Qmask);
% from this step, I need to gather blobs (regridded as well as in physical
% space) for the four quadrants of interest.
% Each quadrant structure needs to contain the following information:
%
for i = 1:NQ
    QInfo(i).cloudInfo = get_subset(cloudInfo_sub, 3, Qmask{i});
    QInfo(i).windInfo = get_subset(windInfo_sub, 3, Qmask{i});
    QInfo(i).SSTInfo = get_subset(SSTInfo_sub, 3, Qmask{i});
    QInfo(i).blobprops = get_subset(blobprops_sub, 2, Qmask{i});

     % cutouts:
     for iv = 1:length(COvarNames)
         VN = COvarNames{iv};
         QInfo(i).cutouts.(VN) = cutouts_sub.(VN)(Qmask{i});
     end
     
     % distance:
     [QInfo(i).dist_metric, Qsid] = sort(r_norm_sel(Qmask{i}),'descend');
     CFa_norm_tmp = CFa_norm_sel(Qmask{i});
     QInfo(i).CFa_norm = CFa_norm_tmp(Qsid);
     
     Conv_norm_tmp = Conv_norm_sel(Qmask{i});
     QInfo(i).Conv_norm = Conv_norm_tmp(Qsid);
     clear CFa_norm_tmp Conv_norm_tmp
     
     blobID_tmp = blobID_sel(Qmask{i});
     QInfo(i).blobID = blobID_tmp(Qsid);
     
     % do the sorting and plotting for each cartesian quadrant.

     % use Qsid to sort the blobInfo.
     QInfo(i).cloudInfo = get_subset(QInfo(i).cloudInfo, 3, Qsid);
     QInfo(i).windInfo = get_subset(QInfo(i).windInfo, 3,Qsid);
     QInfo(i).SSTInfo = get_subset(QInfo(i).SSTInfo, 3,Qsid);
     QInfo(i).blobprops = get_subset(QInfo(i).blobprops, 2, Qsid);
%      QInfo(i).blobprops.mask = QInfo(i).blobprops.mask(Qsid);
%      QInfo(i).blobprops.mask_coord = QInfo(i).blobprops.mask_coord(Qsid);
%      QInfo(i).blobprops.blob_cen = QInfo(i).blobprops.blob_cen(:, Qsid);
%      QInfo(i).blobprops.blob_EqvRadDgr = QInfo(i).blobprops.blob_EqvRadDgr(Qsid);
     
     
     for iv = 1:length(COvarNames)
         VN = COvarNames{iv};
         QInfo(i).cutouts.(VN) = QInfo(i).cutouts.(VN)(Qsid);
     end
end

%% ready to plot the results:
if 1==0
for ifig = 1:4
    
    
    figsvdirsub = [figsvdir filesep 'examine_indiv_features'];
    if ~exist(figsvdirsub, 'dir')
        mkdir(figsvdirsub)
    end
    
        hfig = plot_blob_stamps(QInfo(ifig), ifig, coord, figsvdirsub);    % this will be a nested function


    %figname = ['Q' num2str(ifig) '_indivBlob_SST_CldFrac_WindConv.jpg'];
    %xc_savefig(gcf, figsvdir, figname, [0 0 12 12]);
    
end
end

varargout{1} = selected_frac;
varargout{2} = QInfo;
end

%% uitlity functions:
function hfig = plot_blob_stamps(QInfo, ifig, coord, figsvdir)

%% 2. plotting the results:
ncase = length(QInfo.dist_metric);

% default;
nrow = 6; ncol = 6;

varlist = {'SSTa_cutouts','cloudfrac_highfreq','winddiv_highfreq'};
var_stdname_list = {'SST anom', 'cloud fraction anom', 'wind divergence'};
category_transformed = {'SSTInfo', 'cloudInfo','windInfo'};
category_physical = {'SSTa','CFa','winddiv'};


%vrange = {[-0.5, 0.5],[-0.005,0.005],[-0.5, 0.5],[-0.015, 0.015]}; % SST, cldfrac, div, cldfrac_transformedGrid, div_transformedGrid
vrange = {[-0.5, 0.5],[-0.10, 0.10],[-0.01,0.01]};



figcnt =0;
pos = customize_subplot_size(nrow, ncol, 0.02, 0.05);
for ic =1: ncase
    if mod(ic, nrow)==1
        figcnt = figcnt+1;
        hfig = figure(ifig); clf;
    end
    
    irow = mod(ic,nrow);
    if irow==0
        irow = nrow;
    end
    
    
    %% in the regridded coordinate
    % cloud anomaly relative to climatology &  % wind divergence:
    for icol = 1:3
        VN = varlist{icol};
        catname = category_transformed{icol};
        
        ip = (irow-1)*ncol + icol;
        val = QInfo.(catname).(VN)(:,:,ic);
        
        hsub(ip) = subplot(nrow, ncol, ip);
        pcolor(coord.XX, coord.YY, val);shading flat;
        colorbar;
        colormap(redblue);
        caxis(vrange{icol});
        hold on;
        circle(0,0,1,'k',2);
        var_stdname = var_stdname_list{icol};
        
        
        titlestr = {[var_stdname];'in normalized space'};
        
        if icol==1
            maxSSTa = max(val(:));
            titlestr{3} = ['SSTa_{max} = ' num2str(maxSSTa,'%3.2f')];
        elseif icol ==2
            titlestr{3} = ['NormMetric = ' num2str(QInfo.CFa_norm(ic), '%4.2f')];
        elseif icol==3
            titlestr{3} = ['NormMetric = ' num2str(QInfo.Conv_norm(ic), '%4.2f')];
        end
        
        
        
        if irow ==1
            title(titlestr); 
        else
            title(titlestr([3]));
        end
           
        axis('square');
        set(hsub(ip), 'pos', pos{ip});

    end

    
    %% in original physical coordinate (need to have the associated cloudInfo, windInfo ready for the case..)
    for icol = 4:6
        ip = (irow-1)*ncol + icol;
        
        hsub(ip) =subplot(nrow, ncol, ip);
        % try
        iv = mod(icol,3);
        if iv ==0
            iv =3;
        end
        
        % get the right information from QInfo:
        CN = category_physical{iv};
        tmp = QInfo.cutouts.(CN)(ic);
        
        VN = varlist{iv};
        valvec = tmp.(VN)(:);
        
        if contains(CN, 'SST')
            lonvec = tmp.LON_cutouts(:);
            latvec = tmp.LAT_cutouts(:);
        else
            lonvec = tmp.lon(:);
            latvec = tmp.lat(:);
        end
        
        scatter(lonvec, latvec,30,  valvec,'filled','marker','s');
        colorbar;
        colormap(redblue);
        caxis(vrange{iv});
        hold on;
        
        % plot the SST feature contours here:
        contour(QInfo.blobprops.mask_coord(ic).lon, QInfo.blobprops.mask_coord(ic).lat, ...
            QInfo.blobprops.mask{ic}, 'k');
        xb = QInfo.blobprops.blob_cen(1,ic);
        yb = QInfo.blobprops.blob_cen(2,ic);
        rb = QInfo.blobprops.blob_EqvRadDgr(ic);
        axis([xb-3*rb, xb+3*rb, yb-3*rb, yb+3*rb]);
        
        % plot wind direction in an arrow:
        %mean_wdir = mean(atan2(windsel(ic).vwnd, windsel(ic).uwnd),'all','omitnan');
        %mean_uwnd = mean(QInfo.cutouts.winddiv(ic).uwnd(:),'omitnan');
        %mean_vwnd = mean(QInfo.cutouts.winddiv(ic).vwnd(:), 'omitnan');
        mean_uwnd = QInfo.blobprops.ave_wspd(ic) * cosd(QInfo.blobprops.ave_wdir(ic));
        mean_vwnd = QInfo.blobprops.ave_wspd(ic) * sind(QInfo.blobprops.ave_wdir(ic));
        
        
        hold on
        quiver(xb, yb, mean_uwnd*0.15, mean_vwnd*0.15, 'k', 'linestyle','--','linewidth',1.2);
        var_stdname = var_stdname_list{iv};
        if irow ==1
            title({[var_stdname];'in phsyical space'});
        end
        circle(xb,yb, rb,'c',2);
        axis('square');
        
        set(hsub(ip), 'pos', pos{ip});
        %end
    end
    
    
    %pause

    % save figure to designated location:
    if mod(ic,nrow)==0 || ic ==ncase
        stid = num2str(ic-nrow+1, '%3.3d'); edid = num2str(ic, '%3.3d');
        figname = ['Quadrant' num2str(ifig, '%2.2d') '_C' stid '-C' edid '_SST_cloudfrac_windiv_Fig' num2str(figcnt,'%2.2d') '.jpg'];
        disp(figname);
        xc_savefig(gcf, figsvdir, figname, [0 0 13.5 12]);
    end
end

end


