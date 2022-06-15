function [sigmask, p2D_sub, samplesz] = rank_sum_field_significance_test_with_FDR(blobInfo, varname, subrange,  alpha_FDR, TestTailSide, visualize )
% Purpose: this function is used to test hypothesis that two samples are
% different in location (e.g., mean) using rank sum nonparametric testing. 
% the field significance is further evaluated by controling the false
% dicovery rate (FDR) as discussed in Wilks' statistic text book along with
% an article. 
%
% Null hypothesis: on average, warm and cold blobs have the same effect on
% the variables of interest. 
%
% Test method: one of the classical nonparametric tests for location :
% rank-sum test

if contains(varname, 'wind')
    catname = 'windInfo';
elseif contains(varname, 'cloud')
    catname = 'cloudInfo';
else
    catname = 'SSTInfo';
end

xwarm = blobInfo.warm.(catname).(varname);
xcold = blobInfo.cold.(catname).(varname);

% prepare data above for the rank-sum test:
% remove features where NaN is everywhere:
x1 = xwarm;
x2 = xcold;
exids = [];
for i = 1:size(xwarm,3)
    nanmask = isnan(xwarm(:,:,i));
    if all(nanmask)
        exids = [exids,i];
    end
end
x1(:,:,exids) = [];

exids = [];
for i = 1:size(xcold,3)
    nanmask = isnan(xcold(:,:,i));
    if all(nanmask)
        exids = [exids,i];
    end
end
x2(:,:,exids) = [];

% now pool the two samples together treating them as the same (the null
% hypothesis) and rank them.
n1 = size(x1,3); n2 = size(x2,3);
ntot = n1 + n2;

xnull = cat(3, x1,x2);
% at each grid point, sort the 911 data sample and obtain the rank
[~,sid] = sort(xnull,3);
%tag_sorted = tag(sid);
[~, rnk] = sort(sid,3);   % this describe the rankings of the original data xnull, so it follows already the order of the original data;
 
% rearrange size of arrays;
xnull = permute(xnull, [3,1,2]);
xnull_2D = reshape(xnull, ntot, []);
rnk = permute(rnk, [3,1,2]);
r_2D = reshape(rnk, ntot,[]);

% check to see if there are duplicate values, if so, adjust the rank (rnk);
cnt = 0;
for i = 1:size(xnull_2D,2)
    xnull_slice = xnull_2D(:,i);
    [unique_data, uid] = unique(xnull_slice,'stable'); % output not sorted
    dupids = setdiff(1:ntot, uid);  % duplicate ids;
    if ~isempty(dupids)
        % record how many duplication
        cnt = cnt +1;
        num_dup(cnt) = length(dupids);
    end
    for j = 1:length(dupids)
        locids = find(xnull_slice == xnull_slice(dupids(j)));
        rmean = mean(r_2D(locids,i));
        r_2D(locids,i) = rmean;
    end
end

% now compute R1 and R2 based on the warm and cold blob labels;
R1 = sum(r_2D(1:n1,:),1);
R2 = sum(r_2D(n1+1:n1+n2,:),1);

% get the U statistics:
U1 = R1 - n1/2*(n1+1);
% now find the pvalue base on the null hypothesis distribution
% (the null distribution of the Mann-Whitney U-statistics approximately Guassian);
% only 7967 locations out of 121*121 locations have duplication, and the
% maximum duplication instance is 7 7/911 0.0077 (quite small, use Eq 5.24b
% in Wilks' statistic book)
mu_U = (n1*n2)/2;
sigma_U = sqrt(n1*n2*(n1+n2+1)/12);   % this is standard divation.

Guassian_z = (U1 - mu_U)./sigma_U;    % this is a standard normal distribution, mean = 0 and stdv = 1;

% decide whether or not to do a 1-side or 2-sided P values..
% get a gussian distribution pdf or cdf to find p values to determine
% significance.

% 
mu_z = 0;
sigma_z = 1;
%zx = [-4.5:0.01:4.5];
maxGz = ceil(max(abs(Guassian_z)));
zx = [-maxGz:0.01:maxGz];
cdf_zdist = normcdf(zx);    % this will by default create the standard normal distribution.
zcdf_obsv = interp1(zx, cdf_zdist, Guassian_z);

% use two sided test to get the p values:
if strcmp(TestTailSide, 'single')  % P(Z<=z): if x1 is expected to be smaller than x2;
    p = zcdf_obsv;          % take left tail.
    posmask = Guassian_z>0;  % take right tail (Z>=zobsv)
    % take Z<=z
    p(posmask) = 1-zcdf_obsv(posmask);
%elseif strcmp(TestTailSide,'right')
   % p = ();
elseif strcmp(TestTailSide,'both')
    zcdf_obsv = interp1(zx, cdf_zdist, abs(Guassian_z));
    p = 2*(1-zcdf_obsv);
end


p2D = reshape(p, size(xwarm,1), size(xwarm,2));

p2D_sub = p2D(subrange.y, subrange.x);
samplesz.warm = n1;
samplesz.cold = n2;

% call the FDR method.
sigmask = multiple_hypothesis_test_with_FDR(p2D_sub, alpha_FDR);

if visualize
    [coordsub.XX, coordsub.YY] = meshgrid(1:size(p2D_sub,2), 1:size(p2D_sub,1));
    figure(1); clf;
    pcolor(coordsub.XX, coordsub.YY, mean(x1(subrange.y, subrange.x,:),3,'omitnan')); shading flat; colorbar;
    hold on
    scatter(coordsub.XX(sigmask), coordsub.YY(sigmask), 20, 'x', 'MarkerFaceAlpha',1,'MarkerEdgeColor','c');
    sigmask2 = p2D_sub<0.05;
    scatter(coordsub.XX(sigmask2), coordsub.YY(sigmask2), 10, 'o', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor','k');
    circle(0,0,1,'k',2);
    %caxis([-0.04 0.04]);
    %caxis([-5e-4, 5e-4])
    %axis([-1.5, 1.5, -1.5, 1.5])
    colormap(redblue);
    hold off;

end
return
