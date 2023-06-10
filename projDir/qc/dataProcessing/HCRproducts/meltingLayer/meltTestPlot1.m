%% Plot 1

newInds=1:round(length(dataShort.time)/2000):length(dataShort.time);

% Resample for plotting
newDBZ=dataShort.DBZ_MASKED(:,newInds);
newLDR=dataShort.LDR(:,newInds);
velPlot=dataShort.VEL_MASKED;
velPlot(dataShort.elevation>0)=-velPlot(dataShort.elevation>0);
newVEL=velPlot(:,newInds);
newVELdiff=velDiff(:,newInds);
newASL=dataShort.asl(:,newInds);
newTime=dataShort.time(newInds);

newProb=meltProb(:,newInds);
newProb(newProb<0.1)=nan;

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible',showPlot);

ax1=subplot(4,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
view(2);
colMapDBZ(sub1);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
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
xlim([dataShort.time(1),dataShort.time(end)]);
title('LDR')
grid on
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
xlim([dataShort.time(1),dataShort.time(end)]);
title('VEL')
grid on
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
xlim([dataShort.time(1),dataShort.time(end)]);
title('VEL diff')
grid on
ax4.Position=[0.06 0.05 0.87 0.21];

linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer1_',datestr(dataShort.time(1),formatOut),'_to_',datestr(dataShort.time(end),formatOut)],'-dpng','-r0');
