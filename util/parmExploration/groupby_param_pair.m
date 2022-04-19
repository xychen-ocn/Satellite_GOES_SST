function [gbIDs, pedges] = groupby_param_pair(dataIn, subset_ids, param_pair)
% Purpose: group data by input parameter(s);
% Reference: (xarray groupby function)
%  - visualize composite?? 
%

% define parameters:
p1 = param_pair{1}; p2 = param_pair{2};
%p1_edges = pair_edges{1}; p2_edges = pair_edges{2};


% get parameter values from all the features:
for i = 1:2
    PN = param_pair{i};
    if ~strcmp(PN, {'mean_wspd', 'max_SSTgrad'})
         tmp = [dataIn.properties];
         val_tmp = [tmp.(PN)];
         pval.(PN) = val_tmp(subset_ids);
    else
        if strcmp(PN, 'mean_wspd')
            tmp = [dataIn.cloudInfo];
             val_tmp = [tmp.mean_wspd];
             pval.(PN) = val_tmp(subset_ids);
        else
            % find max_SSTgrad...
            tmp = [dataIn.SSTInfo];
            for j = 1:length(tmp)
                tmp(j).max_SSTgrad = max(abs(tmp(j).eSSTgrad(:)));
            end
            val_tmp = [tmp.max_SSTgrad];
            pval.(PN) = val_tmp(subset_ids);
            
        end
    end
    
    % get the edges for 4 qartiles;
    pedges.(PN) = prctile(pval.(PN), [0:25:100]);

    
end

% now what: start grouping via percentile:
clear gids sstids
p1_gtag=discretize(pval.(p1), pedges.(p1));
for id = 1:length(pedges.(p1))-1
    % loop through different bins of wind speed
    ids_lev1{id} = find(p1_gtag==id);
    
    % for each p1 bin, sort blobs further by p2:
    %SSTidx = discretize(gravel_SSTmax(gids{id}), SSTa_edges);
    p2_gtag = discretize(pval.(p2)(ids_lev1{id}), pedges.(p2));
    
    for k = min(p2_gtag):max(p2_gtag)
        tmp = find(p2_gtag==k);
        if ~isempty(tmp)
            gbIDs{id, k} = tmp;
        end
    end
    %
end

    




return