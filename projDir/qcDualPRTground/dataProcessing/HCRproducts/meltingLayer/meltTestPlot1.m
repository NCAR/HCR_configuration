%% Plot 1

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible',showPlot);

ax1=subplot(4,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
view(2);
colMapDBZ(sub1);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('Reflectivity and melting layer')
grid on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.87 0.21];
ax1.SortMethod='childorder';

% LDR

ax2=subplot(4,1,2);
hold on;
surf(newTime,newASL./1000,newLDR,'edgecolor','none');
view(2);
clim([-25 -5]);
ax2.Colormap=jet;
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('LDR')
grid on
box on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.525 0.87 0.21];

% VEL

ax3=subplot(4,1,3);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ax3.Colormap=velCols;
caxis([-8 8]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('VEL')
grid on
box on
set(gca,'xticklabel',[])
ax3.Position=[0.06 0.287 0.87 0.21];

ax4=subplot(4,1,4);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ax4.Colormap=velCols;
caxis([-8 8]);
scatter(timeMat(maskForPlot==1),aslForMask(maskForPlot==1)./1000,2.5,'filled','MarkerFaceColor','k');
%scatter(data.time,data.iceLev./1000,4,'filled','MarkerFaceColor','g');
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('Melting layer')
grid on
box on
ax4.Position=[0.06 0.05 0.87 0.21];
ax4.SortMethod='childorder';

linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer1_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
