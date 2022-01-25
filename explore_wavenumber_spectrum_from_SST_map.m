% test computing 2-dimensional Fourier spectrum from the Satellite SST
% field with missing data;
clear all; close all;
% add path to the function I have to compute 2D spectrum
addpath('/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_Cloud_aggregation_analysis/bin');

figsvdir = './Figs';
if ~exist(figsvdir,'dir')
    mkdir(figsvdir)
end

% load in satellite GOES data:
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
files = dir([datadir filesep 'GOES_SST*.mat']);
filenames = {files.name};
cld_thres_array = [0.2:0.2:0.8];
cnt2 = 1;

clear Sk ave_Sk ave_Sk_wgt
for j = 1:length(filenames)
    matFN = filenames{j};
    tmp = strsplit(matFN,'_');
    dateID = tmp{end}(1:5);
    
    
    load([datadir filesep matFN]);
    
    disp(['--> working on ' dateID]);
    
    % deal with missing data: fill with the mean value;
   % cnt=1;
    clear Sk_SST LOgive_SST Sk_cloud LOgive_cloud cloudy_pixel_pcen
    intp_method = {'mean','linear'};
    
    for it = 1:24
        % take subset again to exclude islands...
        xstart = -59;
        xmask = GOES_ATOMIC.lon(:,it)>xstart;
        datasub.lon = double(GOES_ATOMIC.lon(xmask,it));
        datasub.lat = double(GOES_ATOMIC.lat(:,it));
        
        [LON, LAT]= meshgrid(datasub.lon, datasub.lat);
        
        SST_map = GOES_ATOMIC.sea_surface_temperature(xmask,:,it)';
        figure(1); clf;
        subplot(3,1,1)
        pcolor(datasub.lon, datasub.lat, SST_map);
        shading flat;
        title('raw data');
        
        mean_SST0 = mean(SST_map,[1,2],'omitnan');
        cloud_flag = isnan(SST_map);
        cloud_pcen = numel(find(cloud_flag==1))/numel(cloud_flag);
        cloud_pcen_all(it+24*(j-1))=cloud_pcen;
        
        
        cloudy_pixel_pcen=[];
        for k = 1:length(cld_thres_array)
            
            cld_thres = cld_thres_array(k);
            cnt=1;
            
            if ~isnan(mean_SST0) && cloud_pcen< cld_thres
                disp(['   --> testing cloud fraction thres' num2str(cld_thres) ]);
                cloudy_pixel_pcen(cnt) = cloud_pcen;
                
                SST_map_tmp = SST_map;
                for im = 1:length(intp_method)
                    MN =  intp_method{im};
                    
                    if strcmp(MN,'mean')
                        SST_map_tmp(cloud_flag==1) = mean_SST0;
                        
                    else                        
                        SST_intp = scatteredInterpolant(LON(~cloud_flag), LAT(~cloud_flag), SST_map(~cloud_flag),intp_method{im});                        
                        % use inpaint option (a different function)
                        SST_map_tmp(cloud_flag==1)=SST_intp(LON(cloud_flag), LAT(cloud_flag));
                    end
                    SST_map_intped.(MN) = SST_map_tmp;
                    
                    figure(1);
                    subplot(3,1,im+1)
                    pcolor(datasub.lon, datasub.lat, SST_map_intped.(MN));
                    shading flat;
                    title(['cloudy pixels filled with ' MN ' interpolated SST']);
                    
                    % compute spectrum
                    if ~strcmp(MN,'mean')
                        mean_SST = mean(SST_map_intped.(MN),[1,2],'omitnan');
                        SST_anom.(MN) = SST_map_intped.(MN) - mean_SST;
                    else
                        SST_anom.(MN) = SST_map_intped.(MN) - mean_SST0;
                    end
                    [wn_bins, Sk(j,k).SST.(MN)(cnt,:), ~] = compute_wavenumber_spectrum_from_2Dimages(datasub,SST_anom.(MN));
                    
                    %pause(0.5)
                    % cloud only
                    [wn_bins, Sk(j,k).cloud(cnt,:), ~] = compute_wavenumber_spectrum_from_2Dimages(datasub,cloud_flag);
                    %                 %pause(0.5)
                    %
                    % do the weighted averages here:
                    
                     
                    ave_Sk(j,k).(MN) = mean(Sk(j,k).SST.(MN), 1);
                    wgt_matrix =repmat(1-cloudy_pixel_pcen',1,size(wn_bins,2));
                    ave_Sk_wgt(j,k).(MN) = sum(wgt_matrix.*Sk(j,k).SST.(MN),1)./sum(1-cloudy_pixel_pcen);
                    
                    
                end
                cnt=cnt+1;
                
            end
        end   %  loop of cloud threshold
        
    end % loop of the SST map
    
    % do a weighted average of the wavenumber spectrum:
    % weighted by the amount of cloudy pixel:
%     clear ave_Sk ave_Sk_wgt
%     if exist('Sk_SST','var')
%         % ave_Sk = (sum(w(i)*Sk(i))/sum(w(i));
%         for im = 1:2
%             MN = intp_method{im};
%             ave_Sk.(MN) = mean(Sk_SST.(MN), 1);
%             ave_Sk_wgt.(MN) = sum(repmat(1-cloudy_pixel_pcen',1,size(wn_bins,2)).*Sk_SST.(MN),1)./sum(1-cloudy_pixel_pcen);
%         end
%         wvlen_bins = 2*pi./wn_bins /1E3;
%         
%         figure(10); clf;
%         %plot(wvlen_bins(1,:), ave_Sk,'linewidth',1.1);
%         hold on;
%         plot(wvlen_bins(1,:), ave_Sk_wgt.linear,'-r','linewidth',1.2);
%         plot(wvlen_bins(1,:), ave_Sk_wgt.natural,'-m','linewidth',1.2);
%         plot(wvlen_bins(1,:), mean(Sk_cloud, 1),'-k','linewidth',1.3);
%         %plot(wvlen_bins(1,:), Sk_SST(1,:),'.-c');
%         %plot(wvlen_bins(1,:), Sk_cloud(1,:),'.-m');
%         y_ideal =50* wn_bins(1,:).^(-5/3);
%         y_ideal2 = 1* wn_bins(1,:).^(-2);
%         y_ideal3 = 1*10^(-4)* wn_bins(1,:).^(-3);
%         
%         hl(1)=plot(wvlen_bins(1,:), y_ideal, ':b');
%         hl(2)=plot(wvlen_bins(1,:), y_ideal2, '--b');
%         hl(3)=plot(wvlen_bins(1,:), y_ideal3, '-.b');
%         
%         
%         set(gca,'xscale','log','yscale','log','xdir','reverse');
%         xlim([10^0, 10^3.5]);
%         %xlabel('wavenumber (rad/m)');
%         xlabel('wavelength (km)');
%         ylabel('S(k) (K^2/km)');
%         grid on
%         set(gca,'fontsize',14)
%         title(dateID)
%         legend('ave SST spec (linear)','ave SST spec (natural)','ave cloud spec','k^{-5/3}','k^{-2}','k^{-3}');
%         hold off
%         xc_savefig(gcf,figsvdir,['GOES16_SST_wnspec_' dateID '_interp_compared.jpg'],[0 0 8 6]);
  %  end
end


% plot 1D spectrum computed with different SST filling method and with
% different cloud fraction:
% remove zeros..
ave_Sk_cloud=nan(26,4,49);
for k=1:4
    for i= 1:26
        if ~isempty(Sk(i,k).cloud)
        ave_Sk_cloud(i,k,:) = mean(Sk(i,k).cloud,1,'omitnan');
        end
    end
end

% average by day
wvlen_bins = 2*pi./wn_bins /1E3;
colorn={'k','r'};
figure;
for k = 1:length(cld_thres_array)
    for im = 1:2
        MN = intp_method{im};
        ave_Sk_alldays =[ave_Sk(:,k).(MN)];
        ave_Sk_tmp = reshape(ave_Sk_alldays, length(wn_bins),[])';  % put along the row direction first.
        
        Sk_ensm.(MN)=mean(ave_Sk_tmp,1);
        
        wave_Sk_alldays =[ave_Sk_wgt(:,k).(MN)];
        wave_Sk_tmp = reshape(wave_Sk_alldays, length(wn_bins),[])';  % put along the row direction first.
        Sk_wensm.(MN) = mean(wave_Sk_tmp,1);
        
        Sk_cloud_ensm = mean(squeeze(ave_Sk_cloud(:,k,:)),1,'omitnan');
        
        subplot(2,2,k)
        plot(wvlen_bins(1,:), Sk_ensm.(MN),'-','linewidth',1.2,'color',colorn{im});
        hold on;
        hl(im)=plot(wvlen_bins(1,:), Sk_wensm.(MN),'-','linewidth',1.2,'marker','*','color',colorn{im});
    end
    
   
    hl(3)=plot(wvlen_bins(1,:), Sk_cloud_ensm,'-c','linewidth',1.1); %,'color',[0.5 0.5 0.5]);
    %plot(wvlen_bins(1,:), Sk_SST(1,:),'.-c');
    %plot(wvlen_bins(1,:), Sk_cloud(1,:),'.-m');
    %y_ideal =50* wn_bins(1,:).^(-5/3);
    %y_ideal2 = 1* wn_bins(1,:).^(-2);
    y_ideal3 = 1*10^(-4)* wn_bins(1,:).^(-3);
    
    %hl(1)=plot(wvlen_bins(1,:), y_ideal, ':b');
    %hl(2)=plot(wvlen_bins(1,:), y_ideal2, '--b');
    hl(4)=plot(wvlen_bins(1,:), y_ideal3, '-.b');
    
    
    set(gca,'xscale','log','yscale','log','xdir','reverse');
    xlim([10^0, 10^3.5]);
    %xlabel('wavenumber (rad/m)');
    xlabel('wavelength (km)');
    ylabel('S(k) (K^2/km)');
    grid on
    set(gca,'fontsize',14)
    title(['cloudy fraction threshold:' num2str(cld_thres_array(k))])
    legend(hl,{'SST spec (mean SST)','SST spec (linear interped)',...
          'cloud spec','k^{-3}'},'fontsize',11);
    hold off
end
xc_savefig(gcf,figsvdir,['GOES16_SST_26day_ensemble_wnspec.jpg'],[0 0 12 10]);

% plot cloudy pixel histogram:
figure;
histogram(cloud_pcen_all,[0:0.05:1], 'Normalization','count');
% hold on;
% histogram(cloud_pcen_cond,[0:0.05:1], 'Normalization','count','FaceAlpha',0.5);

xlabel('Cloudy Pixel Fraction')
ylabel('Count')
title('All SST maps from 26 days');
set(gca,'fontsize',14)
xc_savefig(gcf,figsvdir,['GOES16_SST_CloudyFrac.jpg'],[0 0 8 6]);
