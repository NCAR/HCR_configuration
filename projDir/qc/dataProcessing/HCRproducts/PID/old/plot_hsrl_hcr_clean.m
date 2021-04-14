function plot_hsrl_hcr_clean(data,pid_comb,backscatLog,cscale_comb,units_str_comb,ylimits)
%Plot hsrl and hcr pid
f2=figure('DefaultAxesFontSize',12,'Position',[400 300 1000 1000]);

s1=subplot(4,1,1);
fig1=surf(data.time,data.asl,data.HCR_DBZ,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-30 15]);
colorbar;
ylabel('Altitude (km)');
title('Radar reflectivity (dBZ)');
colormap(s1,jet)

s2=subplot(4,1,2);
fig1=surf(data.time,data.asl,backscatLog,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
colormap(s2,jet);
caxis([-7 -2]);
colorbar;
ylabel('Altitude (km)');
title(['Lidar backscatter log_{10}[{\beta}] m^{-1} Sr^{-1}']);

s3=subplot(4,1,3);
fig1=surf(data.time,data.asl,data.HSRL_Volume_Depolarization,'edgecolor','none');
view(2);
ylim(ylimits);
caxis([0 0.8]);
colormap(s3,jet);
xlim([data.time(1),data.time(end)]);
colorbar;
ylabel('Altitude (km)');
title(['Lidar depolarization ratio']);

s4=subplot(4,1,4);
fig1=surf(data.time,data.asl,pid_comb,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0.5 9.5]);
colormap(s4,cscale_comb);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels=units_str_comb;

ylabel('Altitude (km)');
title(['Particle ID Combined']);
end

