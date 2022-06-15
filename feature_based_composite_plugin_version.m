% Purpose of this script: do composite analysis;
% Usage: called by the main script, data passed in as global variables.

global dataroot casesvdir figsvdir datasvdir
global indiv_blobs_collection downwind_FCN_regridded cutouts_physicalspace
global sample_random_blobs random_blobs rb_downwind_FCN_regridded  


% ============ Start the scripts ============= %

NT = length(indiv_blobs_collection);
FeatureProps = indiv_blobs_collection;     % blob/feature property information.

% update the script here;
% ========================== Obsolete =================================== %
% for i = 1:length(indiv_blobs_collection)
%     feature_based_ds(i).cloudInfo = cloudiness_downwind_FCN{i};        % cloudiness related parameter in a feature centered coordinate.
%     feature_based_ds(i).SSTInfo = SSTcutout_downwind_FCN{i};           % SST, SST anomaly, and effective SST gradient in a feature centered coordinate.
%     feature_based_ds(i).properties = indiv_blobs_collection(i);  
%     feature_based_ds(i).windInfo = winddiv_downwind_FCN{i};
% end
% ========================== Obsolete =================================== %

%% updated to :
for i = 1:length(indiv_blobs_collection)
    feature_based_ds(i).properties = indiv_blobs_collection(i);
    feature_based_ds(i).properties.ave_windInfo = downwind_FCN_regridded(i).ave_windInfo;
    %feature_based_ds(i).SSTInfo = downwind_FCN_regridded(i).SST;   % tmporary until the error is corrected
    feature_based_ds(i).SSTInfo = downwind_FCN_regridded(i).SSTInfo;
    feature_based_ds(i).cloudInfo = downwind_FCN_regridded(i).cldfrac;
    feature_based_ds(i).windInfo = downwind_FCN_regridded(i).winddiv;
    feature_based_ds(i).cutouts_cloud = cutouts_physicalspace(i).cldfrac;
    feature_based_ds(i).cutouts_wind = cutouts_physicalspace(i).winddiv;
    feature_based_ds(i).cutouts_SST = cutouts_physicalspace(i).SSTInfo;
    
    if sample_random_blobs
    rb_based_ds(i).properties = random_blobs(i);
    rb_based_ds(i).SSTInfo = rb_downwind_FCN_regridded(i).SSTInfo;
    rb_based_ds(i).cloudInfo = rb_downwind_FCN_regridded(i).cldfrac;
    rb_based_ds(i).windInfo = rb_downwind_FCN_regridded(i).winddiv;
    end
       
end

% the data above is from a large domain, now I would like to confine my
% serach region to the following:
% total feature number = 6245; (cold + warm);
search_region = [-59, -48; 8, 18];
selIds.all = select_blobs_within_roi(feature_based_ds, search_region);

% do not worry about different atmospheric regimes for now. 
regime_grped_features = group_daily_blobs_by_largescale_U10_and_LTS_conditions(feature_based_ds);

if sample_random_blobs
    selIds_rb.all = select_blobs_within_roi(rb_based_ds, search_region);
    regime_grped_rb = group_daily_blobs_by_largescale_U10_and_LTS_conditions(rb_based_ds);
end

  % get selection indices for each regime group:
regime_names = fieldnames(regime_grped_features);
for i = 1:length(regime_names)
    RN = regime_names{i};
    blobsIn = regime_grped_features.(RN);
    % the input structure has ? layers: 
    %  -1. properties, cloudInfo, and SSTInfo;
    %  -2. data variables in each of the above fields
    selIds.(RN) = select_blobs_within_roi(blobsIn, search_region);
    
    if sample_random_blobs
        rbsIn = regime_grped_rb.(RN);
        selIds_rb.(RN) = select_blobs_within_roi(rbsIn, search_region);
    end
end

%% Figure 1:
%% 0. show blob statistics: 
% boxplot show radius, eccentricity, and background SST, for all period,
% and 4 different regimes.
ATOMIC_feature_charac_2month = get_feature_characteristics_for_selected_times(feature_based_ds, selIds.all);
for i = 1:4
    RN = regime_names{i};
    ATOMIC_feature_charac_regimes(i) = get_feature_characteristics_for_selected_times(regime_grped_features.(RN), selIds.(RN));
end

if sample_random_blobs
    ATOMIC_feature_charac_2month = get_feature_characteristics_for_selected_times(rb_based_ds, selIds_rb.all);
    for i = 1:4
        RN = regime_names{i};
        ATOMIC_feature_charac_regimes(i) = get_feature_characteristics_for_selected_times(regime_grped_rb.(RN), selIds_rb.(RN));
    end
end
%% make the following script a function to make plots (compare blobs stats)
% make 3 matrices to store corresponding characteristic data. 
% each column stores a group (all + 4 regimes)
% add computation for mean value.
blob_types = fieldnames(ATOMIC_feature_charac_2month.EqvRad_km);
pos = customize_subplot_size(2,3,0.1, 0.1);
figure(2);clf;
for irow = 1:2
    BT = blob_types{irow};
    EqvRad = [ATOMIC_feature_charac_2month.EqvRad_km.(BT)'; ...
        ATOMIC_feature_charac_regimes(4).EqvRad_km.(BT)'; ...
        ATOMIC_feature_charac_regimes(1).EqvRad_km.(BT)'; ...
        ATOMIC_feature_charac_regimes(2).EqvRad_km.(BT)'; ...
        ATOMIC_feature_charac_regimes(3).EqvRad_km.(BT)'];
    
    Eccen = [ATOMIC_feature_charac_2month.Eccen.(BT)'; ...
        ATOMIC_feature_charac_regimes(4).Eccen.(BT)'; ...
        ATOMIC_feature_charac_regimes(1).Eccen.(BT)'; ...
        ATOMIC_feature_charac_regimes(2).Eccen.(BT)'; ...
        ATOMIC_feature_charac_regimes(3).Eccen.(BT)'];
    
    
    AxisRatio = [ATOMIC_feature_charac_2month.AxisRatio.(BT)'; ...
        ATOMIC_feature_charac_regimes(4).AxisRatio.(BT)'; ...
        ATOMIC_feature_charac_regimes(1).AxisRatio.(BT)'; ...
        ATOMIC_feature_charac_regimes(2).AxisRatio.(BT)'; ...
        ATOMIC_feature_charac_regimes(3).AxisRatio.(BT)'];

    
    SSTbg = [ATOMIC_feature_charac_2month.SSTbg.(BT)'; ...
        ATOMIC_feature_charac_regimes(4).SSTbg.(BT)'; ...
        ATOMIC_feature_charac_regimes(1).SSTbg.(BT)'; ...
        ATOMIC_feature_charac_regimes(2).SSTbg.(BT)'; ...
        ATOMIC_feature_charac_regimes(3).SSTbg.(BT)'];
    
    % check the mean values:
    mean_EqvRad = [mean(ATOMIC_feature_charac_2month.EqvRad_km.(BT), 'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(4).EqvRad_km.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(1).EqvRad_km.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(2).EqvRad_km.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(3).EqvRad_km.(BT),'omitnan')];
    
    mean_Eccen = [mean(ATOMIC_feature_charac_2month.Eccen.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(4).Eccen.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(1).Eccen.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(2).Eccen.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(3).Eccen.(BT),'omitnan')];
    
    
    mean_AxisRatio = [ mean(ATOMIC_feature_charac_2month.AxisRatio.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(4).AxisRatio.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(1).AxisRatio.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(2).AxisRatio.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(3).AxisRatio.(BT),'omitnan')];

    
    mean_SSTbg = [mean(ATOMIC_feature_charac_2month.SSTbg.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(4).SSTbg.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(1).SSTbg.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(2).SSTbg.(BT),'omitnan'); ...
        mean(ATOMIC_feature_charac_regimes(3).SSTbg.(BT),'omitnan')];
    %%
    
    g1 = repmat({'All'}, size(ATOMIC_feature_charac_2month.EqvRad_km.(BT)'));
    g2 = repmat({'Fish'}, size(ATOMIC_feature_charac_regimes(4).EqvRad_km.(BT)'));
    g3 = repmat({'Gravel'}, size(ATOMIC_feature_charac_regimes(1).EqvRad_km.(BT)'));
    g4 = repmat({'Flowers'}, size(ATOMIC_feature_charac_regimes(2).EqvRad_km.(BT)'));
    g5 = repmat({'Sugar'}, size(ATOMIC_feature_charac_regimes(3).EqvRad_km.(BT)'));
    g = [g1;g2;g3;g4;g5];
    
    
    % fix plot to have the same ylimit.
    hsub1 = subplot(2,3,1+3*(irow-1));
    % feature equivalent radius
    boxplot(EqvRad, g);
    grid on
    hold on
    plot([1:5], mean_EqvRad,'*m');
    ylabel('Eqv. Radius (km)');
    set(gca,'fontsize',11);
    ylim([0, 160]);
    set(gca,'ytick',[0:20:160]);
    set(gca,'yMinorGrid','on','YMinorTick','on');
    set(hsub1,'pos',pos{1+3*(irow-1)});

    
    % feature eccentricity:
    hsub2 =subplot(2,3,2+3*(irow-1));
    boxplot(AxisRatio, g);
    grid on
    hold on
    plot([1:5], mean_AxisRatio,'*m');
    ylabel('Major to Minor Axis Ratio')
    set(gca,'fontsize',11);
    ylim([0, 8]);
    set(hsub2,'pos',pos{2+3*(irow-1)});
    
    % background SST
    hsub3=subplot(2,3,3*irow);
    boxplot(SSTbg-273.15, g);
    grid on
    hold on
    plot([1:5], mean_SSTbg-273.15,'*m');
    ylabel('Background SST (^{\circ}C)');
    ylim([24.5, 28]);
    set(gca,'yMinorGrid','on','YMinorTick','on');
    set(gca,'fontsize',11);
    set(hsub3,'pos',pos{3*irow});
end
figname = 'warm_cold_feature_characteristics_boxplot.jpg';
xc_savefig(gcf, figsvdir, figname, [0 0 11, 7]);


%% Figure 2:
figsvdirsub = [figsvdir filesep 'composite'];
if ~exist(figsvdirsub,'dir')
    mkdir(figsvdirsub)
end

AllFPs = [feature_based_ds.properties];
AllFPs_stats = [AllFPs.stats_selected];
EqvDiam = [AllFPs_stats.EqvDiam]; 
EqvRad_roi = EqvDiam(selIds.all)*AllFPs(1).L4_xres/1e3 *0.5;  %km

large_mask = EqvRad_roi<= 140 & EqvRad_roi>=40;   % radius > 50km; 40km ~ 65th, 140km: ~99th
selIds.all_large = selIds.all(large_mask);

% update the compoiste script below:
coord.XX = downwind_FCN_regridded(1).XX;
coord.YY = downwind_FCN_regridded(1).YY;

[ATOMIC_composites, samplesz, ~, blobInfo, hfig] = composite_data_in_windaligned_coord_warm_and_cold_updated(feature_based_ds, selIds.all_large,'casename','all','coord',coord);
%[ATOMIC_composites, samplesz, coord, blobInfo] = composite_data_in_windaligned_coord_warm_and_cold(feature_based_ds, selIds.all_large,'casename','all','viewflag',false);

% no need to sample random blobs for now;
%% testing null hypothesis: the averaged change of cloudiness over warm and cold blobs are the same. (warm/cold features have no influence on cloudiness.)
% note that SSTa, eSST are significant already.

QOIs = {'cloudfrac_relanom_2mon','winddiv_highfreq'};
rsel = 2;
stid = find(coord.XX(1,:)==-rsel);
edid = find(coord.XX(1,:)==rsel);
subrange.x = stid:edid;
subrange.y = stid:edid;
alpha_FDR = 0.2;
testside = 'single';
%varname = 'cloudfrac_relanom_2mon';
%varname = 'eSSTgrad';
visflag = true;
for iv = 1:length(QOIs)
    varname = QOIs{iv};
    [sigmask.(varname), p2D_sub.(varname), samplesz] = rank_sum_field_significance_test_with_FDR(blobInfo, varname, subrange, alpha_FDR, testside, visflag);
end

hfig = visualize_composite_with_significant_test(coord, ATOMIC_composites, subrange, samplesz, sigmask);
% save the above figures:
figname = 'composite_warm_vs_cold_alldays.jpg';
xc_savefig(hfig, figsvdirsub, figname, [0 0 12 6]);

%% Figure 2bs: repeat the same operation but for different atmospheric regimes:
for tt = 1:4
    
    RN = regime_names{tt};
    dataIn = regime_grped_features.(RN);
    
    AllFPs = [dataIn.properties];
    AllFPs_stats = [AllFPs.stats_selected];
    EqvDiam = [AllFPs_stats.EqvDiam];
    EqvRad_roi = EqvDiam(selIds.(RN))*AllFPs(1).L4_xres/1e3 *0.5;  %km
    
    large_mask = EqvRad_roi<= 140 & EqvRad_roi>=40;   % radius > 50km; 40km ~ 65th, 140km: ~99th
    subids_large = selIds.(RN)(large_mask);
    
    % update the compoiste script below:
    % coord.XX = downwind_FCN_regridded(1).XX;
    % coord.YY = downwind_FCN_regridded(1).YY;
    
    %[ATOMIC_composites, samplesz, ~, blobInfo, hfig] = composite_data_in_windaligned_coord_warm_and_cold_updated(feature_based_ds, selIds.all_large,'casename','all','coord',coord);
    [grped_composites(tt), group_samplesz(tt),~, RgblobInfo] = composite_data_in_windaligned_coord_warm_and_cold_updated(dataIn, subids_large, 'casename',RN, 'coord',coord,'viewflag',false);
    
    % no need to sample random blobs for now;
    %% testing null hypothesis: the averaged change of cloudiness over warm and cold blobs are the same. (warm/cold features have no influence on cloudiness.)
    % note that SSTa, eSST are significant already.
    
    QOIs = {'cloudfrac_relanom_2mon','winddiv_highfreq'};
    rsel = 2;
    stid = find(coord.XX(1,:)==-rsel);
    edid = find(coord.XX(1,:)==rsel);
    subrange.x = stid:edid;
    subrange.y = stid:edid;
    alpha_FDR = 0.2;
    testside = 'single';
    %varname = 'cloudfrac_relanom_2mon';
    %varname = 'eSSTgrad';
    visflag = true;
    for iv = 1:length(QOIs)
        varname = QOIs{iv};
        [sigmask.(varname), p2D_sub.(varname), samplesz] = rank_sum_field_significance_test_with_FDR(RgblobInfo, varname, subrange, alpha_FDR, testside, visflag);
        sigmask2.(varname) = p2D_sub.(varname)<0.05;
       % pause;
    end
    
    hfig1 = visualize_composite_with_significant_test(coord, grped_composites(tt), subrange, samplesz, sigmask);
    figname = ['composite_warm_vs_cold_ranksum_test_FDR_' RN '.jpg'];
    xc_savefig(hfig, figsvdirsub, figname, [0 0 12 6]);

     hfig2 = visualize_composite_with_significant_test(coord, grped_composites(tt), subrange, samplesz, sigmask2);
    % save the above figures:
    figname = ['composite_warm_vs_cold_ranksum_test_naive_' RN '.jpg'];
    xc_savefig(hfig, figsvdirsub, figname, [0 0 12 6]);
    
    
end

%% Figure 3: view case by case?  
% add a function to gather images with the same features together.
% image that has the right pattern of wind divergence:
% cloud enhancement w/ nagative wind div.

% several type:
[casefrac, blobsIn4Qt] = group_and_check_individual_cases(coord, blobInfo.warm, figsvdir);


% acutally, I can further ask myself what environmental characteristics
% these blobs are associated with. % If I can get the Ids of this blobs
% then I can use that Id to take out informations..
% attempt to plot the maximum SSTa and wind speed
dayType = get_dayType_by_U10_and_LTS;    % 4 regimes named by mesoscale organization; contain dates for each regime.

% check if the blob time falls into 4 regimes roughly...
for i = 1:4
    %rtag{i} = zeros(size(flag));
    for r = 1:4
        RN = regime_names{r};
        flag = ismember(floor(blobsIn4Qt(i).blobprops.time), dayType.(RN));        
        %rtag{i}(flag) = r;
        
        % count the occurrence time of each regime:
        rcount(i,r) = sum(double(flag));
    end
    
end

rcount_for_4Qts_prc = rcount./sum(rcount,2);
QtPrc_for_4Regimes = rcount./sum(rcount,1);

figure(5);clf;
subplot(2,2,1);
bar(rcount,'stacked');
hlgd = legend(regime_names);
ylabel('# of blobs');
set(gca,'xticklabel',{'Q1','Q2','Q3','Q4'});
set(gca,'fontsize',14);
grid on
title('Regime Composition');
set(hlgd,'loc','northeastout');

subplot(2,2,3);
bar(rcount','stacked');
hlgd =legend({'Q1','Q2','Q3','Q4'});
ylabel('# of events');
set(gca,'xticklabel',regime_names);
set(gca,'fontsize',14,'XTickLabelRotation',90);
grid on
title('Cloud-Wind Quadrant Composition');
set(hlgd,'loc','northeastout');


subplot(2,2,2);
bar(rcount_for_4Qts_prc*100,'stacked');
hlgd = legend(regime_names);
ylabel('Percentage (%)');
set(gca,'xticklabel',{'Q1','Q2','Q3','Q4'});
set(gca,'fontsize',14); %,'XTickLabelRotation',90);
grid on
title('Regime Composition');
set(hlgd,'loc','northeastout');

subplot(2,2,4);
bar(QtPrc_for_4Regimes'*100,'stacked');
hlgd = legend({'Q1','Q2','Q3','Q4'});
ylabel('Percentage (%)');
set(gca,'xticklabel',regime_names);
set(gca,'fontsize',14,'XTickLabelRotation',90);
grid on
title('Cloud-Wind Quadrant Composition');
set(hlgd,'loc','northeastout');

figname = 'barplot_composition_of_EnvRegime_for_each_CldFracConv_Quadrant.jpg';
xc_savefig(gcf, figsvdir, figname, [ 0 0 12 8]);





if 1==0
%% test:
%% obsolete:
% for random blobs:
AllFPs = [rb_based_ds.properties];
AllFPs_stats = [AllFPs.stats_selected];
EqvDiam = [AllFPs_stats.EqvDiam]; 
EqvRad_roi = EqvDiam(selIds_rb.all)*AllFPs(1).L4_xres/1e3 *0.5;  %km

large_mask = EqvRad_roi<= 140 & EqvRad_roi>=40;   % radius > 50km; 40km ~ 65th, 140km: ~99th
selIds_rb.all_large = selIds_rb.all(large_mask);

[rb_composites, ~, ~, rbInfo,hfig] = composite_data_in_windaligned_coord_warm_and_cold_updated(rb_based_ds, selIds_rb.all_large,'casename','all','coord',coord);

% use bootstrap to estimate a significant level for the case here. 
comb_rb = cat(3, rbInfo.warm.cloudInfo.cloudfrac_relanom_2mon, rbInfo.cold.cloudInfo.cloudfrac_relanom_2mon);
% for i = 1:1000
%     idx = randperm(size(comb_rb,3), 444);
dw = permute(rbInfo.warm.cloudInfo.cloudfrac_relanom_2mon,[3,1,2]);
dc = permute(rbInfo.cold.cloudInfo.cloudfrac_relanom_2mon,[3,1,2]);

%     bstat_2D(:,:,i) = mean(d,1,'omitnan');
% end
nboot = 1000;
[bstat, bsam] = bootstrp(nboot, @(v)mean(v,1,'omitnan'), dw);
bstat_2D = reshape(bstat, nboot, size(coord.XX,1),[]);
siglev95_2D = prctile(bstat_2D, 95, 1);
siglev05_2D = prctile(bstat_2D, 5, 1);

% siglev95_2D = prctile(bstat_2D, 99, 3);
% siglev05_2D = prctile(bstat_2D, 1, 3);

% try masking the results:
insigmask = ATOMIC_composites.warm.cloudfrac_relanom_2mon<=squeeze(siglev95_2D) & ATOMIC_composites.warm.cloudfrac_relanom_2mon>=squeeze(siglev05_2D) ;

[bstat, bsam] = bootstrp(nboot, @(v)mean(v,1,'omitnan'), dc);
bstat_2D = reshape(bstat, nboot, size(coord.XX,1),[]);
siglev95_2D = prctile(bstat_2D, 95, 1);
siglev05_2D = prctile(bstat_2D, 5, 1);
insigmask2 = ATOMIC_composites.cold.cloudfrac_relanom_2mon<=squeeze(siglev95_2D) & ATOMIC_composites.cold.cloudfrac_relanom_2mon>=squeeze(siglev05_2D) ;

figure
subplot(1,2,1);
pcolor(coord.XX, coord.YY, ATOMIC_composites.warm.cloudfrac_relanom_2mon); shading flat;
colormap(redblue); colorbar; caxis([-0.1 0.1]);
hold on;
circle(0,0, 1,'k',2);
scatter(coord.XX(insigmask), coord.YY(insigmask),2,'MarkerEdgecolor', [0.5 0.5 0.5]);
hold off

subplot(1,2,2);
pcolor(coord.XX, coord.YY, ATOMIC_composites.cold.cloudfrac_relanom_2mon); shading flat;
colormap(redblue); colorbar; caxis([-0.1 0.1]);
hold on;
circle(0,0, 1,'k',2);
scatter(coord.XX(insigmask2), coord.YY(insigmask2),2,'MarkerEdgecolor', [0.5 0.5 0.5]);
hold off

% would this result change if I change the signifcant threshold level?

% compute significance:
sigthres = 97.5;  % so that we can access how different the statistic is from the central 95% of the case. 
% significant level for cloud fraction may need to be shifted first by the mean of the feature state ...
[siglevs, bci]=bootstrap_siglevs_and_CI(blobInfo, rbInfo, sigthres, samplesz);



% save composite results: what does this level mean: 95 out of 100 random
% samples, the cloudiness anomalies from a random field below this
% threshold.
% repeat what I have done using composites:
visualize_composite_with_significant_test(coord, ATOMIC_composites, samplesz, siglevs, bci);


end



