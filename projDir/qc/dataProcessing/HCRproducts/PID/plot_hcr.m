function plot_hcr(data,pid_hcr,cscale_hcr,units_str_hcr,ylimits)
%Plot hcr pid
f3=figure('DefaultAxesFontSize',12,'Position',[400 300 1000 1000]);

s1=subplot(5,1,1);
fig1=surf(data.time,data.asl,data.HCR_DBZ,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-30 15]);
colorbar;
ylabel('Altitude (km)');
title('Reflectivity (dBZ)');
colormap(s1,cid_cmap2)

s2=subplot(5,1,2);
fig1=surf(data.time,data.asl,data.HCR_VEL,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-2 2]);
colorbar;
ylabel('Altitude (km)');
title(['Radial velocity (m s^{-1})']);
colormap(s2,jet)

s3=subplot(5,1,3);
fig1=surf(data.time,data.asl,data.HCR_LDR,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-30 -5]);
colorbar;
ylabel('Altitude (km)');
title(['LDR (dB)']);
colormap(s3,cid_cmap)

s4=subplot(5,1,4);
fig1=surf(data.time,data.asl,data.HCR_WIDTH,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0 2]);
colorbar;
ylabel('Altitude (km)');
title(['Spectrum width (m s^{-1})']);
colormap(s4,cid_cmap)

s5=subplot(5,1,5);
fig1=surf(data.time,data.asl,pid_hcr,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0.5 8.5]);
colormap(s5,cscale_hcr);
cb=colorbar;
cb.Ticks=1:8;
cb.TickLabels=units_str_hcr;
ylabel('Altitude (km)');
title(['Particle ID HCR']);
end

