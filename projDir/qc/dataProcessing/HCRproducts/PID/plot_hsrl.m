function plot_hsrl(data,pid_hsrl,backscatLog,cscale_hsrl,units_str_hsrl,ylimits)
%Plot hsrl and hsrl pid
f1=figure('DefaultAxesFontSize',12,'Position',[400 300 1000 1000]);

colormap(jet);

subplot(4,1,1)
fig1=surf(data.time,data.asl,backscatLog,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-7 -2]);
colorbar;
ylabel('Altitude (km)');
title(['Log10 aerosol backscatter coefficient (m^{-1} sr^{-1})']);

subplot(4,1,2)
fig1=surf(data.time,data.asl,data.HSRL_Volume_Depolarization,'edgecolor','none');
view(2);
ylim(ylimits);
caxis([0 0.8]);
xlim([data.time(1),data.time(end)]);
colorbar;
ylabel('Altitude (km)');
title(['Aerosol Linear depolarization ratio'],'interpreter','none');

subplot(4,1,3)
fig1=surf(data.time,data.asl,data.temp,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
%caxis([260 280]);
colorbar;
ylabel('Altitude (km)');
title(['Temperature (K)'],'interpreter','none');

s4=subplot(4,1,4);
fig1=surf(data.time,data.asl,pid_hsrl,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0.5 7.5]);
colormap(s4,cscale_hsrl);
cb=colorbar;
cb.Ticks=1:7;
cb.TickLabels=units_str_hsrl;
ylabel('Altitude (km)');
title(['Particle ID HSRL'],'interpreter','none');
end

