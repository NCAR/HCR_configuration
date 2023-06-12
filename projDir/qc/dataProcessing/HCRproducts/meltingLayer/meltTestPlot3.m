%% Plot 3

fig3=figure('DefaultAxesFontSize',11,'position',[1700,1300,1500,600],'visible','off');

% Medians

ax1=subplot(2,1,1);
hold on;
scatter(dataShort.time(~isnan(velMax)),maxVelAlt./1000,3.5,'filled','MarkerFaceColor','b');
scatter(dataShort.time(~isnan(ldrMax)),maxLdrAlt./1000,3.5,'filled','MarkerFaceColor','g');
plot(dataShort.time,medFilled./1000,'-r','LineWidth',1.5);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Medians')
grid on
box on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.55 0.87 0.39];

% Diff VEL
ax2=subplot(2,1,2);
hold on;
plot(dataShort.time,stdVel,'-m','LineWidth',1.5);
plot(dataShort.time,stdLdr,'-r','LineWidth',1.5);
ylim([0,400]);
ylabel('Std (m)');
xlim([dataShort.time(1),dataShort.time(end)]);
title('Standard devs')
grid on
box on
ax2.Position=[0.06 0.08 0.87 0.39];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'meltLayer3_',datestr(dataShort.time(1),formatOut),'_to_',datestr(dataShort.time(end),formatOut)],'-dpng','-r0');
