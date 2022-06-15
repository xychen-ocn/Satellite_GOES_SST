function structOut = reorganize_struct_to_array(structIn)

% Purpose: reduce nx1 structures containing fields with the same shape
%          to a 1x1 structure.
% 
% switch nargin
%     case 2
%         n=varargin{1};
% end

if length(structIn)>1
    fields = fieldnames(structIn);
    for iv = 1:length(fields)
        FN = fields{iv};
        if ~any(strcmp(FN, {'XX', 'YY'}))
            
            n = ndims(structIn(1).(FN));
            if numel(structIn(1).(FN))==length(structIn(1).(FN)) || n==3
                % matlab vectors
                n = n-1;
            end
            structOut.(FN) = cat(n+1, structIn.(FN));
        end
    end
    
end


return