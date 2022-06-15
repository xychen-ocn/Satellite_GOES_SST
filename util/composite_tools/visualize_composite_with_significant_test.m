function hfig = visualize_composite_with_significant_test(coord, feature_comp, subrange, sample_size, sigmask)
% Purpose: visualize composite results with significant test added upon it.
%           Do this for all the variables;
% significant test: field significance with rank-sum classical
% nonparametric test
%
% Date : Jun 8, 2022
% loop over different blobs and different variables of interest.

blob_type = fieldnames(feature_comp);
varlist_sigtest = fieldnames(sigmask);
%varlist_sigtest = {'cloudfrac', 'cloudfrac_relanom_2mon', 'winddiv_highfreq'};
coordsub.XX = coord.XX(subrange.y, subrange.x);
coordsub.YY = coord.YY(subrange.y, subrange.x);

%%  make plots
%% start plotting the results to make sense of it.
% loop through variables and construct the significant and confidence
% interval mask for each variable. 
% I have a visualize compoiste function, not usable.

% setting up plotting parameters.
cmap_div = redblue(32);
cmap_jet = jet(32);
varlist = {'SSTa_cutouts','eSSTgrad', 'SSTlap', ...      %  'cloudfrac',
    'cloudfrac_relanom_2mon','winddiv_highfreq'}; %,'cloudfreq_std', 'cldfreq_anom_std'};     % why I still don't have cloudfreq_anom, needs to be added.??
varname = {'SST anom.', '\bf{u} \cdot \nablaSST', '\nabla^2SST', ... %'cloudiness',
    {'cloudiness spatial anomalies'; ' relative to 2month mean'}, {'wind divergence'}}; %{'cloudiness anomalies', 'from "climatology"' }%, 'cloudiness stdv.', 'cloudiness anomaly stdv.'};
varunits = {'K', 'K/hr', 'K/100km^2', '', 's^{-1}'};
colormaps = {cmap_div, cmap_div, cmap_div, cmap_div, cmap_div};  %cmap_jet,
caxis_range_default ={[-0.15, 0.15], [-7.5*10^(-2), 7.5*10^-2], [-0.05, 0.05], [-0.05,0.05], [-5e-4, 5e-4]};    % [0.3, 0.5],, [-0.3,0.3]
xticks = {[-0.15:0.05:0.15], [-6e-2:2e-2:6e-2], [-0.05:0.01:0.05], [-0.05:0.01:0.05], [-5e-4:1e-4:5e-4]};
    
nrow = 2;  ncol=length(varlist);
xspace = 0.08; yspace = 0.08;
pos = customize_subplot_size(nrow, ncol, xspace, yspace);

hfig = figure(10); clf;
set(hfig, 'pos',[339, 410, 1280,600]);
%set(hfig, 'Name', casename);
for ir = 1:nrow
    BT = blob_type{ir};
    
    for icol = 1:ncol      % warm and cold features
        VN = varlist{icol};
        val = feature_comp.(BT).(VN)(subrange.y, subrange.x);
        maxval = max(abs(val(:)));
        %caxis_range = [-maxval, maxval];
        
        
        ip = icol+ncol*(ir-1);
        hsub = subplot(nrow, ncol, ip);   % different variables:
        
        pcolor(coordsub.XX, coordsub.YY, val); shading flat;
        hold on
        circle(0,0,1,'k',2);
        
        % add confidence interval and significant testing if they exist;
        % denote the insignificant locations and locations fall outside of
        % the confidence interval.
        % the following can be improved.
        if ismember(VN, varlist_sigtest)% & icol~=3
            % prep mask from confidence testing.
%             CImask = prep_confidenceInterval_masks(val, bci.(VN).(BT));
%             Sigmask = prep_significantTesting_masks(val, siglevs.(VN).(BT));
            mask = sigmask.(VN);
            sigloc.XX = coordsub.XX(mask);
            sigloc.YY = coordsub.YY(mask);
            scatter(sigloc.XX, sigloc.YY, 2, '.', 'MarkerFaceAlpha',0.6,'MarkerEdgeColor',[0.55 0.55 0.55]);
        end

        cmap = colormaps{icol};
        colormap(hsub, cmap);
        hb = colorbar(hsub);
        set(get(hb,'xlabel'),'String', varunits{icol});
        set(hb,'xtick', xticks{icol});
        caxis(caxis_range_default{icol});
        

        if mod(ip,ncol)==1
            titlestr = {varname{icol} ; ...
                ['samples (' num2str(sample_size.(BT)) ')']};
        else
            titlestr = varname{icol};
        end
        
        title(titlestr,'fontsize',10.5);

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

end



% extra function called within this script:
function SigMask = prep_significantTesting_masks(dataIn, siglevs)
% the significant level here has two bounding threshold:
siglb = siglevs(1,:,:);
sigub = siglevs(2,:,:);
SigMask = dataIn>=squeeze(sigub) | dataIn<=squeeze(siglb) ;

end

function CImask = prep_confidenceInterval_masks(dataIn, bci)

CImask = dataIn >= squeeze(bci(1,:,:))& ...
    dataIn <=squeeze(bci(2,:,:));


end