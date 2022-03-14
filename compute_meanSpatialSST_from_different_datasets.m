% This script is used to compute the SST anomalies from different datasets,
% using different definition of the mean field (daily mean, 3-day, 5-day
% running averages)
% the mean SST (spatial map) will be found by using Fourier spectral
% filtering.
%
%
% date: Jan 26, 2022 (XYC draft)
%       Note: I think spectral filtering may not be the best way to remove
%       smaller scale variations. (Perhaps EOF is the way, need to ask
%       about this.)
%       Feb 01, 2022 (updated LPF to use gaussian filter)
% ------------------------------------------------------------------------%
clear all; clc; close all;

dataroot = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
DatasetNames={'g16','gpblend','ostia'};

%% 1. load in data separately for the entire Jan-Feb ATOMIC period.
%dataFN.g16 = 'GOES_SST_JanFeb2022_ATOMIC_region.mat';
dataFN.g16 = 'GOES_SST_*.mat';
dataFN.gpblend = 'gpblend_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N.mat';
dataFN.ostia = 'ostia_L4_daily_JanFeb_ATOMIC_lon62W_TO_42W_lat5N_TO_25N.mat';


for i =1:length(DatasetNames)
    dataName = DatasetNames{i};
    
    if i==1
        % read g16 daily data:
        g16files = dir([dataroot filesep dataName filesep dataFN.g16]);
        g16fns= {g16files.name};
        for j = 1:length(g16fns)
            load([dataroot filesep dataName filesep g16fns{j}]);
            data.g16(j) = GOES_ATOMIC;
        end
        
    else
        tmp = load([dataroot filesep dataName filesep dataFN.(dataName)]);
        
        fieldN = fieldnames(tmp);
        data.(dataName) =tmp.(fieldN{1});
    end
end

%% 2. compute the mean field in different fashion:
movmean_windows = [3, 5, 7];     % this will remove some days in the blob searching, choose to return only average with the full window size. 
ave_windows = [15, 30];     % not the running ones;
cutoff_scale_km = 600;           % scale larger than which the SST variability is kept as the large scale. 
                                 % need to explore a bit how to define the
                                 % cutoff_scale; (800km is the barotropic Radius of
                                 % deformation with a mixed layer depth of
                                 % 80m in this region; 600km is for 40m mixed layer depth.)
subregion = [-59, -45; 8 20];    % avoid fft complication by land during 

cutoff_scale_km = 1000;
calspec = false;
%% g16 
% --> construct temporal mean SST field by averaging available hrly maps using 1, 3 to 5 day
% running average (1 day average is not good enough to use spectrum
% filtering because of missing data)
tday = unique(floor([data.g16.time_num]));
%SST_QCed(data.g16.quality_level~=5)=NaN;     % the sea_surface temperature above is already at the best quality level

% 1a. find the daily averaged SST and perform low-pass filter to obtain the
%     large scale SST gradient with wavenumber > the Rossby deformation
%     radius estimated at this region.
for it = 1:length(tday)
    t = tday(it);
    %     tmask = data.g16(it).time_num>=t & data.g16.time_num<t+1;
    %     disp(num2str(length(find(tmask))));
    %     if length(find(tmask))~=24
    %         disp(datestr(data.g16.time_num(tmask)));
    %         pause
    %     end
    SST_QCed = data.g16(it).sea_surface_temperature;
    
    % 1. create daily averaged map for each day
    G16_daily_ave.SST(:,:,it) = mean(SST_QCed,3,'omitnan');
end
G16_daily_ave.time = tday;
G16_daily_ave.lon = data.g16(it).lon;
G16_daily_ave.lat = data.g16(it).lat;
% save this daily averaged data:
datasvdir = [dataroot filesep 'g16/smaller_domain' filesep 'G16_daily_averaged_SST.mat']; 
save(datasvdir,'G16_daily_ave');


% 1b. then, do the running averaged and apply spatial low-pass filtering to
%     obtain the large scale trend in the SST.
lon = double(data.g16(1).lon(:,1));
lat = double(data.g16(1).lat(:,1));

% subregion = [-60, -45; 8 20];    % avoid fft complication by land 
% 
lonmask = lon>=subregion(1,1) & lon<=subregion(1,2);
latmask = lat>=subregion(2,1) & lat<=subregion(2,2);
% 
% G16_movaved.lon_sub = lon(lonmask); G16_movaved.lat_sub = lat(latmask);
G16_movaved.lon = lon; G16_movaved.lat = lat;

for i = 1:length(movmean_windows)
    wz = movmean_windows(i);
    fieldN = ['wndsz_' num2str(wz) 'd'];
    
    SSTin = G16_daily_ave.SST;
    
    %% compute moving averages:
    G16_movaved.(fieldN).SST = movmean(SSTin, wz ,3, 'omitnan', 'Endpoints','discard');      % operate on the time dimension;
    G16_movaved.(fieldN).time = movmean(G16_daily_ave.time, wz, 'Endpoints','discard');
    
    %% low pass filtering to obtain large scale graident:
    %SST_LPF = zeros(size(G16_movaved.(fieldN).SST(lonmask, latmask,:)));
    SST_LPF = zeros(size(G16_movaved.(fieldN).SST));
    SST_LPF = permute(SST_LPF, [2,1,3]);
    
    if 1==1
        clear Sk LOgive
    for it = 1:length(G16_movaved.(fieldN).time)
        SST0 = G16_movaved.(fieldN).SST(:,:,it)';

        if mod(it,10)==0
            checkflag = true;
        else
            checkflag = false;
        end
        
        if calspec
            %% where do I cutoff the scale? this is a question. check the
            %% wave number spectrum to decide?
            % check wavenumber spectrum of the SST to further justify the cutoff scale:
            SST0_demean = SST0 - mean(SST0(:),'omitnan');
            SST0_demean(isnan(SST0_demean)) = 0;
            
            coord.lon = lon(lonmask);
            coord.lat = lat(latmask);
            [wnbins{1}, Sk(it,:), LOgive(it)]=compute_wavenumber_spectrum_from_2Dimages(coord, SST0_demean(latmask, lonmask));
            
            AveWnSpec(1).(fieldN) = mean(Sk,1);
            AveOgiveLen(1).(fieldN) = mean(LOgive,'omitnan');
        end

        % I thnk I may need to use a large map to do this filtering and
        % then take the results from a subdomain?
        SST_LPF(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);    % the size is different here.
                     
%                  test1 =  estimate_LSG_SST(lon(lonmask), lat(latmask), SST0(latmask, lonmask), ...
%                          'method','spectrum','CutoffScale',600, 'checkflag',checkflag);    % the size is different here.
%                      
%                  test2 = estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale',600, 'checkflag',checkflag);    % the size is different here.
%         
%         
%         if checkflag
%             pause
%         end


    end
    end
    
    G16_movaved.(fieldN).SST_LS = SST_LPF;
    
end
close all



% 1c. do the period averaged without sliding the window of averaging.    (averages)
G16_aved.lon = lon; G16_aved.lat = lat;
%G16_aved.lon_sub = lon(lonmask); G16_aved.lat_sub = lat(latmask);
for i = 1:length(ave_windows)
    wz =ave_windows(i);
    fieldN = ['wndsz_' num2str(wz) 'd'];
    
    SSTin = G16_daily_ave.SST;                                % landmask or bad data;
    nt = size(SSTin, 3);
    time_num = G16_daily_ave.time;
    
    ngroup = ceil(nt/wz);
    time_num_full = nan(wz*ngroup,1);
    for ig = 1:ngroup
        if wz*ig <=nt
            idx = 1+(ig-1)*wz:wz*ig;
        else
            idx = 1+(ig-1)*wz:nt;
        end
        G16_aved.(fieldN).SST(:,:,ig) = mean(SSTin(:,:,idx), 3, 'omitnan');
        
    end
    time_num_full(1:nt) = time_num;
    G16_aved.(fieldN).time =reshape(time_num_full, wz, []);
    
    % estimate large scale gradient:
    % large scale filtering:
    %SST_LPF = zeros(size(G16_aved.(fieldN).SST(lonmask, latmask,:)));
    SST_LPF = zeros(size(G16_aved.(fieldN).SST));

    SST_LPF = permute(SST_LPF, [2,1,3]);         % lon-lat-time to lat-lon-time for matlab.
    clear Sk LOgive
    for it = 1:ngroup
        SST0 = G16_aved.(fieldN).SST(:,:,it)';
        lon = double(data.g16(1).lon(:,1));
        lat = double(data.g16(1).lat(:,1));
 
        if calspec
            %% where do I cutoff the scale? this is a question. check the
            %% wave number spectrum to decide?
            SST0_demean = SST0 - mean(SST0(:),'omitnan');
            SST0_demean(isnan(SST0_demean)) = 0;
            
            coord.lon = lon(lonmask);
            coord.lat = lat(latmask);
            [wnbins{1}, Sk(it,:), LOgive(it)]=compute_wavenumber_spectrum_from_2Dimages(coord, SST0_demean(latmask, lonmask));
            AveWnSpec(1).(fieldN) = mean(Sk,1);
            AveOgiveLen(1).(fieldN) = mean(LOgive,'omitnan');
        end
        
        SST_LPF(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',false);
        
    end
    
    G16_aved.(fieldN).SST_LS = SST_LPF;           % low pass filted large scale, XX-day mean SST 
    
end


%% gplend & ostia: first find the mean SST_tot, and then apply the spatial filter.
% --> 3-day running averaged, spatially low-passed filtered large scale SST field.
% --> 5-day running averaged, spatially low-passed filtered large scale SST field

for k = 2:3
    product = DatasetNames{k};
    
    lon = double(data.(product).lon(:,1));
    lat = double(data.(product).lat(:,1));
   
    lonmask = lon>=subregion(1,1) & lon<=subregion(1,2);
    latmask = lat>=subregion(2,1) & lat<=subregion(2,2);

    data_movaved.lon = lon; data_movaved.lat = lat;
    %data_movaved.lon_sub = lon(lonmask); data_movaved.lat_sub = lat(latmask);
    
    % moving average:
    for i = 1:length(movmean_windows)
        wz = movmean_windows(i);
        fieldN = ['wndsz_' num2str(wz) 'd'];
        
        SSTin = data.(product).analysed_sst;
        mask = data.(product).mask;
        SST_masked = SSTin;
        SST_masked(mask==0) = NaN;                                             % landmask or bad data;

        % compute moving averages:
        data_movaved.(fieldN).SST = movmean(SST_masked, wz ,3, 'omitnan', 'Endpoints','discard');      % operate on the time dimension;
        data_movaved.(fieldN).time = movmean(data.(product).time_num, wz, 'Endpoints','discard');
        
        
        % low pass filtering to obtain large scale graident:
        %SST_LPF = zeros(size(data_movaved.(fieldN).SST(lonmask, latmask,:)));
        SST_LPF = zeros(size(data_movaved.(fieldN).SST));
        SST_LPF = permute(SST_LPF, [2,1,3]);
        
        if 1==1
            clear Sk LOgive
        for it = 1:length(data_movaved.(fieldN).time)
            SST0 = data_movaved.(fieldN).SST(:,:,it)';
            
            if mod(it,10)==0
                checkflag = true;
            else
                checkflag = false;
            end
            
            if calspec
                % check wavenumber spectrum of the SST to further justify the cutoff scale:
                SST0_demean = SST0 - mean(SST0(:),'omitnan');
                SST0_demean(isnan(SST0_demean)) = 0;
                
                coord.lon = lon(lonmask);
                coord.lat = lat(latmask);
                [wnbins{k}, Sk(it,:), LOgive(it)]=compute_wavenumber_spectrum_from_2Dimages(coord, SST0_demean(latmask, lonmask));
                AveWnSpec(k).(fieldN) = mean(Sk,1);
                AveOgiveLen(k).(fieldN) = mean(LOgive,'omitnan');
            end
            
            %% where do I cutoff the scale? this is a question. check the
            %% wave number spectrum to decide?
           % cutoff_scale_km = 1000;
            SST_LPF(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);
            
        end
        end
        
        data_movaved.(fieldN).SST_LS = SST_LPF;
        
    end
    eval([upper(product) '_movaved = data_movaved;']);
    clear data_movaved

    
    % 15-day and monthly averaged (fixed)
    data_aved.lon = lon; data_aved.lat = lat;
    %data_aved.lon_sub = lon(lonmask); data_aved.lat_sub = lat(latmask);
  
    for i = 1:length(ave_windows)
        wz =ave_windows(i);
        fieldN = ['wndsz_' num2str(wz) 'd'];
        
        SSTin = data.(product).analysed_sst;
        mask = data.(product).mask;
        SST_masked = SSTin;
        SST_masked(mask==0) = NaN;                                             % landmask or bad data;
        nt = size(SSTin, 3);
        time_num = data.(product).time_num;
        
        ngroup = ceil(nt/wz);
        time_num_full = nan(wz*ngroup,1);
        %% temporal averaging:
        for ig = 1:ngroup
            if wz*ig <=nt
                idx = 1+(ig-1)*wz:wz*ig;
            else
                idx = 1+(ig-1)*wz:nt;
            end
            data_aved.(fieldN).SST(:,:,ig) = mean(SST_masked(:,:,idx), 3);

        end
        time_num_full(1:nt) = time_num;
        data_aved.(fieldN).time =reshape(time_num_full, wz, []);
        
        % estimate large scale gradient:
        %% large scale filtering:
        %SST_LPF = zeros(size(data_aved.(fieldN).SST(lonmask, latmask,:)));
        SST_LPF = zeros(size(data_aved.(fieldN).SST));
        SST_LPF = permute(SST_LPF, [2,1,3]);
        
        if 1==1
            clear Sk LOgive
        for it = 1:ngroup
            SST0 = data_aved.(fieldN).SST(:,:,it)';
            lon = double(data.(product).lon(:,1));
            lat = double(data.(product).lat(:,1));
            checkflag = false;
            
            if calspec
                % check wavenumber spectrum of the SST to further justify the cutoff scale:
                SST0_demean = SST0 - mean(SST0(:),'omitnan');
                SST0_demean(isnan(SST0_demean)) = 0;
                
                coord.lon = lon(lonmask);
                coord.lat = lat(latmask);
                [wnbins, Sk(it,:), LOgive(it)]=compute_wavenumber_spectrum_from_2Dimages(coord, SST0_demean(latmask, lonmask));
                AveWnSpec(k).(fieldN) = mean(Sk,1);
                AveOgiveLen(k).(fieldN) = mean(LOgive,'omitnan');
            end
            
            %% where do I cutoff the scale? this is a question. check the
            %% wave number spectrum to decide?
            SST_LPF(:,:,it) =estimate_LSG_SST(lon, lat, SST0, 'method','spectrum','CutoffScale', cutoff_scale_km, 'checkflag',checkflag);
            
            
        end
        end
        
        data_aved.(fieldN).SST_LS = SST_LPF;

    end

    eval([upper(product) '_aved = data_aved;']);
    clear data_aved

    
end


%% check the averaged wavenumber spectrum for different products: to determine if the 800km cutoff scale is reasonable.
% 600km seems reasonable.
% make a plot to show that 600km is reasonable, but 500 is also okay -->
% test sensitivity later on.
if 1==0

    
    colors{1} = winter(5);
    colors{2} = spring(5);
    colors{3} = turbo(5);
    
    y_ideal =50* wnbins{1}.^(-5/3);
    y_ideal2 = 0.5*1* wnbins{1}.^(-2);
    y_ideal3 = 2*10^(-4)* wnbins{1}.^(-3);
    
    
    figure(10);clf;
    
    for i = 1:3
        fieldN = fieldnames(AveWnSpec(i));
        wvlen_bins = 2*pi./wnbins{i} /1E3;
        
        for j = 1:length(fieldN)
            FN = fieldN{j};
            Sk = AveWnSpec(i).(FN);
            plot(wvlen_bins, Sk,'color', colors{i}(j,:),'linewidth',1.2,'linestyle','-');
            hold on;
        end
    end
    %hl(1)=plot(wvlen_bins, y_ideal, ':b');
    hl(2)=plot(2*pi./wnbins{1} /1E3, y_ideal2, '--b');
    hl(3)=plot(2*pi./wnbins{1} /1E3, y_ideal3, '-.b');
    plot(2*pi./wnbins{1} /1E3, y_ideal3*0.02, '-.b');
    plot([400, 400],[10^5, 10^13],'--k');
    plot([600, 600],[10^5, 10^13],'--k');
    plot([800, 800], [10^5, 10^13], '--k');
    set(gca,'xscale','log','yscale','log','xdir','reverse');
    xlim([10^0, 10^3.5]);
    ylim([10^5, 10^13]);
    xlabel('wavelength (km)');
    ylabel('S(k) (units?)');
    grid on
    lgd = legend(hl(2:3),{'k^{-2}','k^{-3}'});
    set(lgd, 'location','best');
    set(gca,'fontsize',14);
    xc_savefig(gcf,'Figs','SSTspectrum_from_3products_with_different_temporal averaging.jpg',[0 0 10 8]);
end

note = {'dimensions:';'SST_LSF (nlat x nlon x ntime), SST (nlon x nlat x ntime)'};
%% 3. save the mean field for other script to use.
% save it to each dataset's own directory;
datasvFN = ['TemporalAveraged_LowPassFiltered_thres' num2str(cutoff_scale_km) 'km_SST.mat'];
save([dataroot '/g16/' datasvFN],'G16_movaved','G16_aved','note','-v7.3');
save([dataroot '/gpblend/' datasvFN],'GPBLEND_movaved', 'GPBLEND_aved','note');
save([dataroot '/ostia/' datasvFN],'OSTIA_movaved', 'OSTIA_aved','note');

if calspec
    save([dataroot filesep 'AveSpectrum_of_SSTmap_from_differentsource.mat'], 'AveWnSpec','AveOgiveLen','wnbins');
end

