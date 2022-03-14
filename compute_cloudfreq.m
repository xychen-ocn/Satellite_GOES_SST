function CF= compute_cloudfreq(SSTmaps)

dims =size(SSTmaps);
n = length(dims);

nt = size(SSTmaps,n);
cldflag = isnan(SSTmaps);
CF = sum(cldflag,n)/nt;



end