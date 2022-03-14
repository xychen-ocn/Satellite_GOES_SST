function [wn_bins, Sk, LOgive] = compute_wavenumber_spectrum_from_2Dimages_v2(datasub, cldmask, varargin)

% Input: datasub --> a matlab structure
%        cldmask --> a logical matrix. (1= cloud, 0=not cloud)
% 
switch nargin 
    case 2
        units = 'latlon';
    case 3
        units = varargin{1};
end
if strcmp(units, 'latlon')
    midlat = 0.5*(datasub.lat(1)+datasub.lat(end));
    Lx = (datasub.lon - datasub.lon(1))*111E3*cosd(midlat);
    Ly = (datasub.lat - datasub.lat(1))*111E3;

elseif strcmp(units, 'cartesian')
    Lx = datasub.x;
    Ly = datasub.y;
    
else
    disp('invalid units for the input grid');
    return
    
end
dx = diff(Lx);
dy = diff(Ly);

%% from here on, the code can be built into a matlab function. 

%% perform fft2:
Y = fft2(cldmask);
S = abs(fftshift(Y)).^2;

% 1. wavenumber spectrum;
[Ny, Nx] = size(cldmask);

% -- 1.a construct wavenumber space:
if mod(Nx,2)==0
  ix = linspace(-Nx/2,Nx/2-1,Nx);
else
  ix = linspace(-floor(Nx/2),floor(Nx/2),Nx);
end

if mod(Ny,2)==0
  iy = linspace(-Ny/2,Ny/2-1,Ny);
else
  iy = linspace(-floor(Ny/2),floor(Ny/2),Ny);
end

kx = 2*pi/(Nx*mean(dx,'omitnan')) .* ix;   % the last kx is the nyquist frequency
ky = 2*pi/(Ny*mean(dy,'omitnan')) .* iy; 
[KX, KY] = meshgrid(kx,ky);
Kmag = sqrt(KX.^2 + KY.^2);


% construct a wavenumber space that is a square (nx = ny) using the
% largest of (dx and dy),
dcom = min(mean(dx,'omitnan'),mean(dy,'omitnan'));   %
N = min(Nx, Ny);
if N ==Nx
    ind = ix;
else
    ind = iy;
end
kcom = 2*pi/(N*dcom) .* ind;
wn_bins_e = sqrt(kcom.^2*2);
wn_bins = 0.5*(wn_bins_e(1:end-1)+ wn_bins_e(2:end));

% [KXC, KYC]= meshgrid(kcom);
% SC = interp2(KX, KY, S, KXC, KYC);
Sk= zeros(length(wn_bins),1);
for i = 1:length(wn_bins)
    mask = sqrt(KX.^2 + KY.^2)>=wn_bins_e(i) & sqrt(KX.^2 + KY.^2)<wn_bins_e(i+1); 
    Sk(i) = sum(S(mask),'omitnan');
end
%wn_bins = Kmag_mean;

Stot=[];
for i=2:length(wn_bins)
    Stot(i) = trapz(wn_bins(1:i), Sk(1:i));
end

%  --- 3.a find the critical wavenumber using Ogive. --- 
% [uqid,ia, ic] = unique(log10(Stot));
% kc_exponent = interp1(Stot(ia), log10(wn_bins(ia)), 1/3*Stot(end));
% kc = 10^kc_exponent;   % 2/3 of the variance comes from wavenumber larger than kc.
% LOgive = 2*pi/kc/1E3;
LOgive=NaN;
% 

%% make plots to show the spectrum of the input data:
figure(10); clf;
% show the input cloudmask:
hsub(1)=subplot(2,2,1);
[XX, YY] = meshgrid(Lx,Ly);
pcolor(XX, YY, double(cldmask));shading flat;
colormap(hsub(1), gray)
colorbar;

% show the spectrum in kx-ky space:
hsub(2) = subplot(2,2,2);
pcolor(KX, KY, S);shading flat;
ylim([-0.25,0.25]*10^(-3))
xlim([-0.25,0.25]*10^(-3))
xlabel('kx (rad/m)');
ylabel('ky (rad/m)');
colormap(hsub(2), parula)
colorbar;
caxis([0 10^6])



% show the omnidirectional spectrum in the wavelength space:
subplot(2,2,[3,4]);
y_ideal =50* wn_bins.^(-5/3);
y_ideal2 = 1* wn_bins.^(-2);
y_ideal3 = 2*10^(-4)* wn_bins.^(-3);
wvlen_bins = 2*pi./wn_bins /1E3;

plot(wvlen_bins, Sk,'-k');
hold on;
hl(1)=plot(wvlen_bins, y_ideal, ':b');
hl(2)=plot(wvlen_bins, y_ideal2, '--b');
hl(3)=plot(wvlen_bins, y_ideal3, '-.b');

yrange = get(gca,'ylim');
%hl(4)=plot(2*pi./[kc,kc]/1E3,[10^5, 10^10],'--r');

set(gca,'xscale','log','yscale','log','xdir','reverse');
%set(gca,'xtick',fliplr([10^3, round(2*pi/kc/1E3), 10^2, 10^1, 10^0]));
%text(2*pi/kc/1E3, 10^5*5, num2str(round(2*pi./kc/1E3)));
xlim([10^0, 10^3.5]);
%xlabel('wavenumber (rad/m)');
xlabel('wavelength (km)');
ylabel('S(k) (units?)');
grid on
lgd = legend(hl,{'k^{-5/3}','k^{-2}','k^{-3}'}); %,'k_c (Ogive)'});
%set(lgd, 'location','best');


return 