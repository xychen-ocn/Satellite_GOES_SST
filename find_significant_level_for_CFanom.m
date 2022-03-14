clear all; clc; %close all;

datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
svdataFN = 'ATOMIC_JanFeb_SSTmaps.mat';
load([datadir filesep svdataFN]);
load([datadir filesep 'cloudiness_freq_ATOMIC_JanFeb_updated.mat'],'CF_3dmovmean_2monthsAve','CF_JanFeb');

ROI = [-59, -48; 8, 20];
[LON, LAT] = meshgrid(lon,lat);
lonmask = lon>=ROI(1,1) & lon<=ROI(1,2);
latmask = lat>=ROI(2,1) & lat<=ROI(2,2);

mask = lonmask&latmask;

SSTsub_all = SST_all(latmask, lonmask, :);
clear SST_all;


CF_ref = CF_3dmovmean_2monthsAve(latmask, lonmask);                  % this is actually 5d ave, but the differences are small.
sample_info.sampleSize=[405, 209, 143, 617];      % gravel, flowers, sugar, fish
% start with a simpler independent sampling for each regime.
N0 = size(SSTsub_all,3);
Ntest = 1000;
for ig = 2:4     % four regimes:
    for it =1:Ntest
        ridx = randperm(N0, sample_info.sampleSize(ig));
        SSTmaps=SSTsub_all(:,:,ridx);
        %times_here = time_all_operate(ridx);
        % compute cloudiness frequency:
        tmp = compute_cloudfreq_and_anom(SSTmaps, CF_ref);
        indv_random_CFanom{ig}(:,:,it) = tmp.CF_anom;
        clear tmp
        if mod(round(it/Ntest*100), 10)==0
            disp(['Group' num2str(ig) ' finished ' num2str(round(it/Ntest*100)) '% sampling...'])
        end
    end
    
    % compute variance:
    indv_CFa_variance(ig,:) = var(indv_random_CFanom{ig}, 0, [1,2]);                  % length:Ntest
    indvtest_var_thres(ig,:) = prctile(indv_CFa_variance(ig,:), [90, 95]);

end
testsvdir = './sigTest';
if ~exist(testsvdir)
    mkdir(testsvdir);
end

for ig = 1:4
    
    % compute variance:
    tmp = prctile(indv_random_CFanom{ig}, [90,95,99],[1,2]);
    indv_CFa_extremes(ig,:,:) = squeeze(tmp);    % ngroup xnprc x Ntest
    
    % find percentile from distribution
    indvtest_extreme_thres(ig,:) = prctile(indv_CFa_extremes(ig,:,:), 95,3);
end

note{1} = 'var_thres: 90 and 95th prctile';
note{2} = 'extreme_thres: 95th prctile, where extreme are defined as the 90, 95, 99 prctile of in the CFanom samples';
save([testsvdir filesep 'significant_threshold_from_individual_sampling_in_4regimes.mat'],'indv_CFa_variance','indvtest_var_thres','indvtest_extreme_thres')
save([testsvdir filesep 'significant_threshold_from_individual_sampling_in_4regimes.mat'],'note','-append');

% these values are quite small.
figure
for ig = 1:4
    subplot(2,2,ig)
    histogram(indv_CFa_variance(ig,:),20);
end

%% set up to do the random grouping:
draw_order_combinations = {[1:4], [1,2,4,3], [1,3,2,4], [1,3,4,2], [1,4,2,3], [1,4,3,2], ...
[2,1,3,4], [2,1,4,3], [2,3,1,4], [2,3,4,1], [2,4,1,3], [2,4,3,1], ...
[3,1,2,4], [3,1,4,2], [3,2,1,4], [3,2,4,1], [3,4,1,2], [3,4,2,1], ...
[4:-1:1],  [4,3,1,2], [4,2,1,3], [4,2,3,1], [4,1,2,3], [4,1,3,2]};

for it = 1:Ntest
    SSTsub_all_operate = SSTsub_all;
   % time_all_operate = time_all;
 
   j = randi(24);
   draw_order = draw_order_combinations{j};
   
    for id = 1:4
        ig = draw_order(id);
        if id==1
            N = N0;
        end
        ridx = randperm(N, sample_info.sampleSize(ig));
        SST_maps=SSTsub_all_operate(:,:,ridx);
        SSTsub_all_operate(:,:,ridx)=[];       % remove the sample already drew.
        %time_all_operate(ridx)=[];
        N = size(SSTsub_all_operate,3);
        
        % compute cloudiness frequency:
        cloud_flag = isnan(SST_maps);
        CF = sum(double(cloud_flag),3)./size(cloud_flag,3);
        CFa = (CF  - CF_ref)./CF_ref;

        random_grouping_CFanom{ig}(:,:,it) = CFa;
        
    end
    if mod((it/Ntest*100), 10)==0
            disp(['Group' num2str(ig) ' finished ' num2str(round(it/Ntest*100)) '% sampling...'])
    end
end


%% find the variance threshold:
for ig = 1:4
    
    % compute variance:
    %grp_CFa_variance(ig,:) = var(random_grouping_CFanom{ig}, 0, [1,2]);
    tmp = prctile(random_grouping_CFanom{ig}, [90,95,99],[1,2]);  % ngroup x nprc x Ntest
    grp_CFa_extremes(ig,:,:) = squeeze(tmp);
    
    % find percentile from distribution
    grptest_var_thres(ig,:) = prctile(grp_CFa_variance(ig,:), 95);      % use 95 as the threshold.
    grptest_extreme_thres(ig,:) = prctile(grp_CFa_extremes(ig,:,:), 95,3);
end

figure
for ig = 1:4
    subplot(2,2,ig)
    histogram(grp_CFa_variance(ig,:),20);
end
save([testsvdir filesep 'significant_threshold_from_random_regrouping_in_4regimes.mat'],'grp_CFa_variance','grptest_var_thres','grptest_extreme_thres')
save([testsvdir filesep 'significant_threshold_from_random_regrouping_in_4regimes.mat'],'note','-append');


% note: manually checked that the variance of CF_anom from the four regimes
% all above the 95th var threshould determined from random sampling.