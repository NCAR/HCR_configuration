%% Plot 2

fig2=figure('DefaultAxesFontSize',11,'position',[1700,1300,1500,1200],'visible',showPlot);

% Diff VEL
ax1=subplot(3,1,1);
hold on;
surf(newTime,newASL./1000,newVELdiff,'edgecolor','none');
view(2);
colDiff1=winter(4);
colDiff2=spring(20);
colDiff=cat(1,flipud(colDiff1),colDiff2);
ax1.Colormap=colDiff;
caxis([0 0.6]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('VEL diff')
grid on
box on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.87 0.21];

% Diff DBZ
ax2=subplot(4,1,2);
hold on;
surf(newTime,newASL./1000,newDBZdiff,'edgecolor','none');
view(2);
ax2.Colormap=jet;
caxis([-1 4]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('DBZ diff')
grid on
box on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.525 0.87 0.21];

% Melt probability
ax3=subplot(4,1,3);
hold on;
surf(newTime,newASL./1000,newProb,'edgecolor','none');
view(2);
caxis([0 1]);
colProb1=winter(round(thresholds.meltProbLow*10*2));
colProb2=spring(round((thresholds.meltProbHigh-thresholds.meltProbLow)*10*2));
colProb3=copper(round((10-thresholds.meltProbHigh*10)*2));
colProb=cat(1,flipud(colProb1),colProb2,colProb3);
ax3.Colormap=colProb;
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('MeltProb')
grid on
box on
set(gca,'xticklabel',[])
ax3.Position=[0.06 0.287 0.87 0.21];

% Melting layer
meltPlot=nan(size(newMeltLayer));
meltPlot(newMeltLayer==9)=3;
meltPlot(newMeltLayer==11)=2;
meltPlot(newMeltLayer==19)=1;
meltPlot(newMeltLayer==21)=0;

ax4=subplot(4,1,4);
hold on;
surf(newTime,newASL./1000,meltPlot,'edgecolor','none');
view(2);
ax4.Colormap=[0,1,1;0.5,0.5,0.5;0,0,0;1,0,1];
caxis([-0.5,3.5]);
cb=colorbar;
cb.Ticks=[0,1,2,3];
cb.TickLabels={'Cold','Melting cold','Melting warm','Warm'};
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('Melting layer')
grid on
box on
ax4.Position=[0.06 0.05 0.87 0.21];

linkaxes([ax1,ax2,ax3,ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer2_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
