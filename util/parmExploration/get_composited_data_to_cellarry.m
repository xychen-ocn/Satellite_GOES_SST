function [comp_cell_warm, comp_cell_cold] = get_composited_data_to_cellarry(dsIn, parm)
% get structre data into cell array for plotting purpose;
%
blob_types = {'warm', 'cold'};
tmp = [dsIn.(parm)];

parm_struc = reshape(tmp, 4,4)';     % nrow x ncol

for ir = 1:4
    for ic = 1:4
        comp_cell_warm{ir, ic} = parm_struc(ir, ic).warm;
        comp_cell_cold{ir, ic} = parm_struc(ir, ic).cold;
    end
end

return
