% multiple_hypothesis_test by controling false discovery rate:
function sigmask = multiple_hypothesis_test_with_FDR(p_field, alpha_FDR)
% Reference: Wilks: "The Stippling shows statistically significant
% gridponits": how research results are routinely overstated and
% over-interpreted, and what to do about it.
% 
%
% when there are strong spatial correlation (like coherent structure), 
% alpha_FDR = 2* alpha_global
%
% Inputs: p_field (2D matrix the same size of the field of a quantity
%         alpha_FDR: false discovery rate level (0~1, 0.10 an arbitrary
%         default)
% Outputs:

Npts = numel(p_field);         % total number of locations

% 1. sort the p values in the ascending order
pvec = reshape(p_field, 1,[]);  % a row vector;
pvec_ascend = sort(pvec, 'ascend');

% 2. define p*FDR: 
%    sometimes there is ps lower than the pFDR at larger Is..
%    the following function needs improvment...
idx = 1:Npts;
eps = 1e-5;
selmask = pvec_ascend-alpha_FDR.*idx./Npts<=eps;
pvec_sel = pvec_ascend(selmask);

if ~isempty(pvec_sel)
    id_sel = find(selmask==1);
    
    % find the largest index that are consecutive.
    edgeIDs = find(diff(id_sel)>=10); % the separation loc;
    if ~isempty(edgeIDs)
        dim1Dist = [edgeIDs(1); diff(edgeIDs)';length(id_sel)- edgeIDs(end)];
    else
        dim1Dist = [length(id_sel)];
    end
    % break this into different segments and select the end index on the
    % longest segment.
    %the edgeIDs store the cutting point of the data record
    if iscolumn(id_sel)
        tmp_cell= mat2cell(pvec_sel,dim1Dist, 1);
    else
        tmp_cell = mat2cell(pvec_sel, 1, dim1Dist);
    end
    
    nseg = length(tmp_cell);
    if nseg>1
        for i = 1:nseg
            len_seg(i) = length(tmp_cell{i});
        end
        [~,mxid] = max(len_seg);
        pseg_sel = tmp_cell{mxid};
        
    else
        pseg_sel = pvec_sel;
        
    end
    
    if ~isempty(pseg_sel)
        
        p_FDR = max(pseg_sel);
        
        % test null hypothesis by comparing p_field and p_FDR
        % sigmask: 1(reject H0, significant)
        sigmask = p_field<=p_FDR;
        
    else
        % none of the p values are smaller than the FDR pvalue. insignificant
        % results everywhere... (H0 null hypothesis can not be rejected
        % everywhere.)
        disp('none of the p values is smaller than pFDR');
        sigmask = false(size(p_field));
    end


else
    disp('none of the p values is smaller than pFDR');
    sigmask = false(size(p_field));
    return
end


    
     

figure(10);clf;
plot(idx, pvec_ascend,'-r');
hold on
plot(idx, alpha_FDR.*idx./Npts,'--k');
if ~isempty(pvec_sel)
    plot(idx(selmask), pvec_ascend(selmask),'*c');
    if ~isempty(pseg_sel)
        plot([0, max(idx(selmask))], [pseg_sel(end), pseg_sel(end)],'--m');
    end
    xlim([0 max(idx(selmask))+10])
    hold off;
end
%
% i=1;
% while i<=Npts
%     crit = pvec_ascend(i)<=alpha_FDR*(i/Npts);
%     if ~crit        
%         pvec_ascend(i:end) = [];
%         break
%     else
%         i = i+1;
%     end
% end



return