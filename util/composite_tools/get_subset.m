function struct_subseted = get_subset(structIn, dimN, subset_ids)
% pupose: subset fields contained in a structure;
fieldn = fieldnames(structIn);
for iv = 1:length(fieldn)
    FN = fieldn{iv};
    if isnumeric(structIn.(FN))
        if length(size(structIn.(FN)))==3
            if length(dimN)==1
                if dimN==3
                    struct_subseted.(FN) = structIn.(FN)(:,:,subset_ids);
                elseif dimN == 2
                    struct_subseted.(FN) = structIn.(FN)(:,subset_ids,:);
                elseif dimN ==1
                    struct_subseted.(FN) = structIn.(FN)(subset_ids,:,:);
                end
                
            elseif length(dimN)==2
                if dimN(1)==1 && dimN(2) ==2
                    dim1_ids = subset_ids{1};
                    dim2_ids = subset_ids{2};
                    struct_subseted.(FN) = structIn.(FN)(dim1_ids,dim2_ids,:);
                end
            end
            
        elseif length(size(structIn.(FN)))==2
            if dimN == 2
                struct_subseted.(FN) = structIn.(FN)(:,subset_ids);
            elseif dimN ==1
                struct_subseted.(FN) = structIn.(FN)(subset_ids,:);
            end
        end
           
        
    elseif iscell(structIn.(FN))
        struct_subseted.(FN) = structIn.(FN)(subset_ids);
%         if isnumeric(subset_ids)
%             for i = 1:length(subset_ids)
%                 struct_subseted.(FN){i} = structIn.(FN){subset_ids(i)};
%             end
%         elseif islogical(subset_ids)
%             subset_locs = find(subset_ids);
%             for i = 1:length(subset_locs)
%                 struct_subseted.(FN){i} = structIn.(FN){subset_locs(i)};
%             end
%         end
        
    elseif isstruct(structIn.(FN))
        if length(structIn.(FN))>1
            struct_subseted.(FN) = structIn.(FN)(subset_ids);
        else
            disp('can not operate on this structural field')
        end
        
    end
        
    
    
end

end