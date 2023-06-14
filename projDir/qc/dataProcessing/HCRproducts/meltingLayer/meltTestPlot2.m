%% Plot 2

fig2=figure('DefaultAxesFontSize',11,'position',[1700,1300,1500,900],'visible',showPlot);

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
ax1.Position=[0.06 0.69 0.87 0.27];

% Diff DBZ
ax2=subplot(3,1,2);
hold on;
surf(newTime,newASL./1000,newDBZdiff,'edgecolor','none');
view(2);
% colDiff1=winter(4);
% colDiff2=spring(20);
% colDiff=cat(1,flipud(colDiff1),colDiff2);
% ax2.Colormap=colDiff;
ax2.Colormap=jet;
caxis([-1 4]);
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('VEL diff')
grid on
box on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.375 0.87 0.27];

% Melt probability
ax3=subplot(3,1,3);
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
ax3.Position=[0.06 0.06 0.87 0.27];

linkaxes([ax1,ax2,ax3],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer2_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
