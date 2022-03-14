ncFN = 'era_ATOMIC_region_averaged_U10_LTS.nc';
ncdisp(ncFN);

U10 = ncread(ncFN,'U10_domainAve');
LTS = ncread(ncFN, 'LTS_domainAve');

figure
plot(LTS, U10,'.k','color',[0.5 0.5 0.5]);
hold on;
plot(LTS_hrly, U10_hrly)
plot([10, 20], [8, 8],'--r','linewidth',1.2,'color');
plot([15, 15], [4, 12],'--r','linewidth',1.2);
title('3 hourly domain averaged data points in Jan, Feb 2020')
set(gca,'fontsize',14);
xlabel('LTS (K)');
ylabel('U_{10} (m/s)');
grid on