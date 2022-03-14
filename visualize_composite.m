function hfigs = visualize_composite(coord, composites, type, samples)


[nrow, ncol] = size(composites);
nrow_max = min(4, nrow);

pos =customize_subplot_size(nrow_max, ncol, 0.05,0.04);
for i = 1:nrow
    for j = 1:ncol
        ip = j+ncol*(i-1);
        ifig = ceil(ip/nrow_max/ncol);
        if ip > nrow_max*ncol
            ip_sub = ip - nrow_max*ncol;
        else
            ip_sub = ip;
        end
        hfigs(ifig) = figure(ifig); 
        if ~isempty(composites{i,j})
            if samples(i,j)>=10
                    hsub = subplot(nrow_max, ncol, ip_sub);

            if strcmp(type, 'domain_mean')
                cmean = mean(composites{i,j},'all','omitnan');
            elseif strcmp(type, 'downwind_mean')
                cmean = mean(composites{i,j}, 1, 'omitnan');
            end
            pcolor(coord.XX, coord.YY, (composites{i,j}-cmean)); shading flat;
            hold on

            circle(0,0,1,'k',2);
            colorbar; colormap(redblue);
            title(['sample (' num2str(samples(i,j)) ')'], 'fontsize',11);
            caxis([-0.1 0.1])
            axis('square');
            
            set(hsub,'pos', pos{ip_sub});
            end
        end
    end
end


return