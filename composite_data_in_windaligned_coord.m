function [composites, coord, samples, hfigs] = composite_data_in_windaligned_coord(dataIn, subset_ids, methodstr)
% purpose: composite the input data (structure) into different subsets
% using the provided indices:

gridX = [-5:0.1:5]; gridY = gridX;
[nrow, ncol] = size(subset_ids);

nrow_max = min(4, nrow);
pos =customize_subplot_size(nrow_max, ncol, 0.03,0.04);

h = findobj('type','figure');

fignum = max([h.Number]);
%nfig= ceil(nrow/4);

for i = 1:nrow
    for j = 1:ncol
        % get subsets out using indices:
        keys = subset_ids{i,j};
        samples(i,j) = length(keys);    % number of composite samples

        if ~isempty(keys)
            subdata = dataIn(keys);
            gridded_data = grid_data_into_windaligned_coord(subdata, gridX, gridY);
            
            if strcmp(methodstr,'aveCF')
                composites{i,j} = mean(gridded_data.cloudyfreq, 3,'omitnan');
            elseif strcmp(methodstr, 'fromCC')
                composites{i,j} = sum(gridded_data.cloudycnt, 3)./sum(gridded_data.samplesz);
            end
%             %% put results:
            ip = j+ncol*(i-1);
            ifig = ceil(ip/nrow_max/ncol);
            if ip > nrow_max*ncol
                ip_sub = ip - nrow_max*ncol;
            else
                ip_sub = ip;
            end
            hfigs(ifig)=figure(ifig+fignum);
            hsub = subplot(nrow_max, ncol, ip_sub);
            pcolor(gridded_data.XX, gridded_data.YY, composites{i,j}); shading flat;
            hold on
            circle(0,0,1,'k',2);
            colorbar; colormap(turbo);
            title(['samples (' num2str(samples(i,j)) ')'], 'fontsize',11);
            caxis([0.1, 0.6])
            axis('square');
            if mod(ip_sub, ncol)==1
            ylabel({'downwind'; 'normalized distance'});
            end
            
            if ceil(ip_sub/ncol)==nrow_max
                xlabel('crxwind normalized distance');
            end
            
            set(hsub,'pos', pos{ip_sub});
        end
           

    end
end

coord.XX = gridded_data.XX;
coord.YY = gridded_data.YY;
       

return