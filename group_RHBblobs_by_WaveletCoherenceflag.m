function grouped_blobs = group_RHBblobs_by_WaveletCoherenceflag(timerange, blobsIn)

blobs_time = [blobsIn.time];

% get blobs catched within the time range:
ndays = size(timerange,2);

idx_sel = [];
for i = 1:ndays
    cond_time = (blobs_time-timerange(1,i)).*(blobs_time-timerange(2,i));
    if any(cond_time<=0)
        idx = find(cond_time<=0);
        % catched blobs within timerange on this day:
        idx_sel = [idx_sel , idx];
    end
end
grouped_blobs=blobsIn(idx_sel);
return
        