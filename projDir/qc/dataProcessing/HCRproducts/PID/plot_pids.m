function plot_pids(data,pid_comb,pid_comb2,cscale_comb,units_str_comb,ylimits)
f4=figure('DefaultAxesFontSize',12,'Position',[400 300 1000 800]);

s1=subplot(2,1,1);
fig1=surf(data.time,data.asl,pid_comb,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([.5 9.5]);
colormap(s1,cscale_comb);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels=units_str_comb;
ylabel('Altitude (km)');
title(['Particle ID Combined']);

s2=subplot(2,1,2);
fig1=surf(data.time,data.asl,pid_comb2,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([.5 9.5]);
colormap(s2,cscale_comb);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels=units_str_comb;
ylabel('Altitude (km)');
title(['Particle ID Combined Direct']);
end