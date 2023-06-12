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
caxis([-25 -5]);
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
caxis([-6 6]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('VEL')
grid on
box on
set(gca,'xticklabel',[])
ax3.Position=[0.06 0.287 0.87 0.21];

linkaxes([ax1 ax2 ax2 ax3],'xy');

% Diff VEL
ax4=subplot(4,1,4);
hold on;
surf(newTime,newASL./1000,newVELdiff,'edgecolor','none');
view(2);
colDiff1=winter(4);
colDiff2=spring(20);
colDiff=cat(1,flipud(colDiff1),colDiff2);
ax4.Colormap=colDiff;
caxis([0 0.6]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('VEL diff')
grid on
box on
ax4.Position=[0.06 0.05 0.87 0.21];

linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer1_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
