% This script will plot the innate likelihood map (indicate location of the
% location of warm anomalies) and the changes in cloudiness frequency
% together.
datadir = '/Users/xchen/Documents/MATLAB/shallow_convection/Obsv/Satellite_GOES_SST';
CF_file = 'cloudiness_freq_ATOMIC.mat';     % can be updated to the entire DJF
Warmblob_likelihood_FN = 'likelihood_of_warmblob.mat';

load([datadir filesep CF_file]);

product = 'g16';                  % 'gpblend'
TempAveTag = '5d-movingMean';     % '3d-movingMean'
SpatailFilterTag ='LPF600km'; % 'LPF600km';    % 'no_spatial_filtering'
%cloudmaskTag = 'no_cloudmask';
cloudmaskTag = '';
if isempty(cloudmaskTag)
    casedir = strjoin({product, TempAveTag, SpatailFilterTag},'_');
else
    casedir = strjoin({product, TempAveTag, SpatailFilterTag, cloudmaskTag},'_');
end

load([datadir filesep 'blob_data' filesep casdir filesep Warmblob_likelihood_FN]);

% setup colors to plot the probability:
levs = [0:0.1:0.8];
for i = 1:length(levs)-1
    pdf_levs{i} = levs(i:i+1);
end

pdf_color0 = flipud(hot(64));
cindx = floor(linspace(5,60, length(pdf_levs)));
pdf_color = pdf_color0(cindx,:);

parula_all = gray(64);
parula_mid = parula_all(16-10:48+10,:);

% plot:
% how to adjust the hue of the colormap;
figure(1);
for tt = 1 %:4
    CloudType = cloud_types{tt};
    nt = sample_info.sampleSize(tt);
    prc = sample_info.samplePrc(tt)*100;
    
    %     if tt~=3
    %         cindx = [1, 28, 40];
    %         pdf_color = pdf_color0(cindx,:);
    %
    %     else
    %         cindx = [1, 40, 62];
    %         pdf_color = pdf_color0(cindx,:);
    %
    %     end
    %
    
    val = day.(CloudType).CF_anom;
    figure(tt);clf;
    % subplot(2,2,tt);
    pcolor(LON, LAT,val); shading flat;
    %colormap((parula_mask));
    %colormap(rdbusub(floor(linspace(20,60,10)),:));
    colormap(parula_mid);
    hb=colorbar;
    %hb.Ticks=[0.1:0.1:0.8];
    caxis([-0.5, 0.5]);
    set(get(hb,'xlabel'),'string','fractional change','fontsize',12);
    %set(hb,'Ticks',[0:0.1:1]);
    hold on;
    [C,h]=contour(LON,LAT,val, [-0.4:0.2:0.4],'color',[0.25 0.25 0.25]);
    %clabel(C,h, [0:0.2:0.4],'labelspacing',800,'color','k','fontsize',12);
    
    for i = 1:length(pdf_levs)
        c = pdf_color(i,:);
        % plot with the conditioned pdf:
        [Cpdf, hpdf] = contour(XBinCen, YBinCen, P_warmblobs_innate.all,pdf_levs{i},'color',c,'linewidth',1.1, 'linestyle','--');
        % now, plot it with the total pdf;
        %[Cpdf, hpdf] = contour(pdf_data.XBinCen, pdf_data.YBinCen, pdf_data.pdf_all',pdf_levs{i},'color',c,'linewidth',1.1, 'linestyle','--');
        
        clabel(Cpdf,hpdf,'color',c,'labelSpacing',400);
        hold on;
    end
    %     [C3,h3]= contour(LON, LAT, SST_LS_trend_fit, [298.5:0.5:300.5],'w','linestyle',':','linewidth',1.8);
    %     clabel(C3,h3,[299:0.5:300],'color','m', 'labelspacing',800,'fontsize',12);
    xlabel('Longitude (^{\circ}E)');
    ylabel('Latitude (^{\circ}N)');
    hold off;
    set(gca,'fontsize',14);
    %title({'Frequency of missing data (due to cloud overhead)'; 'between Jan09 and Feb12 2020'},'fontsize',16);
    title({'Relative change in cloudiness frequency'; [CloudType '-favored condition in daytime ']},'fontsize',16);
    axis('equal');
    xlim([-59, -48]);
    ylim([10, 18]);
    %xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_'
    %CloudType '_v2.jpg'],[0 0 10 8]);  with SST:
    % xc_savefig(gcf,'Figs',['cloudiness_frequency_GOES_SSTretrival_JanFeb2020_' CloudType '_with_warmblob_likelihood_v0.jpg'],[0 0 10 8]);
    %xc_savefig(gcf,'Figs',[CloudType '_cloudiness_frequency_GOES_SSTretrival_JanFeb2020_with_warmblob_pdf_alltime.jpg'],[0 0 10 8]);
    
end
xc_savefig(gcf,figsvdir,['gravel_daytime_relativeChangeCF_with_bloblikelihood_from_' casedir '_grayhot_corrected.jpg'],[0 0 10 8]);