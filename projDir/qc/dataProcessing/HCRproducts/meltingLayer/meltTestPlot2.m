%% Plot 2

fig2=figure('DefaultAxesFontSize',11,'position',[1700,1300,1500,1200],'visible',showPlot);

timeMat=repmat(dataShort.time,size(data.DBZ_MASKED,1),1);

ax1=subplot(4,1,1);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
view(2);
ax1.Colormap=velCols;
caxis([-6 6]);
scatter(timeMat(meltMask==1),dataShort.asl(meltMask==1)./1000,2.5,'filled','MarkerFaceColor','k');
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Mask')
grid on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.87 0.21];
ax1.SortMethod='childorder';

% Melt probability

ax2=subplot(4,1,2);
hold on;
surf(newTime,newASL./1000,newProb,'edgecolor','none');
view(2);
caxis([0 1]);
colProb1=winter(round(meltProbThreshLow*10*2));
colProb2=spring(round((meltProbThreshHigh-meltProbThreshLow)*10*2));
colProb3=copper(round((10-meltProbThreshHigh*10)*2));
colProb=cat(1,flipud(colProb1),colProb2,colProb3);
ax2.Colormap=colProb;
colorbar
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('MeltProb')
grid on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.525 0.87 0.21];

% Medians

ax3=subplot(4,1,3);
hold on;
scatter(dataShort.time(~isnan(velMax)),maxVelAlt./1000,3.5,'filled','MarkerFaceColor','b');
scatter(dataShort.time(~isnan(ldrMax)),maxLdrAlt./1000,3.5,'filled','MarkerFaceColor','g');
plot(dataShort.time,medFilled./1000,'-r','LineWidth',1.5);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Medians')
grid on
set(gca,'xticklabel',[])
ax3.Position=[0.06 0.287 0.87 0.21];

% Diff VEL
ax4=subplot(4,1,4);
hold on;
plot(dataShort.time,stdVel,'-m','LineWidth',1.5);
plot(dataShort.time,stdLdr,'-r','LineWidth',1.5);
ylim([0,400]);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Standard devs')
grid on
ax4.Position=[0.06 0.05 0.87 0.21];

%linkaxes([ax1 ax2 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer2_',datestr(dataShort.time(1),formatOut),'_to_',datestr(dataShort.time(end),formatOut)],'-dpng','-r0');
