%% Plot 2

fig2=figure('DefaultAxesFontSize',11,'position',[1700,1300,1500,600],'visible',showPlot);

timeMat=repmat(timeForMask,size(data.DBZ_MASKED,1),1);

ax1=subplot(2,1,1);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ax1.Colormap=velCols;
caxis([-6 6]);
scatter(timeMat(maskForPlot==1),aslForMask(maskForPlot==1)./1000,2.5,'filled','MarkerFaceColor','k');
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('Mask')
grid on
box on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.55 0.87 0.39];
ax1.SortMethod='childorder';

% Melt probability

ax2=subplot(2,1,2);
hold on;
surf(newTime,newASL./1000,newProb,'edgecolor','none');
view(2);
caxis([0 1]);
colProb1=winter(round(thresholds.meltProbLow*10*2));
colProb2=spring(round((thresholds.meltProbHigh-thresholds.meltProbLow)*10*2));
colProb3=copper(round((10-thresholds.meltProbHigh*10)*2));
colProb=cat(1,flipud(colProb1),colProb2,colProb3);
ax2.Colormap=colProb;
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([newTime(1),newTime(end)]);
title('MeltProb')
grid on
box on
ax2.Position=[0.06 0.08 0.87 0.39];

linkaxes([ax1,ax2],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer2_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0');
