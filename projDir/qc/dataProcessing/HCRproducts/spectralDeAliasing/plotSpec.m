function plotSpec(data,sampleNum,startInd,powerSpecLarge,ylimUpper,powerSpecFilt,showPlot,figdir,saveWaterfall)
f1 = figure('Position',[100 500 1500 1100],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

subplot(1,2,1)
surf(1:sampleNum*5,data.range./1000,powerSpecLarge,'EdgeColor','none');
view(2);

xlim([1,sampleNum*5]);
ylim([0,data.range(end)./1000]);

ylim([0 ylimUpper])
xticks(1:sampleNum:sampleNum*5);

xlabel('Sample number')
ylabel('Range (km)')
title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))

caxis([0 20])
colorbar

subplot(1,2,2)
surf(1:sampleNum*5,data.range./1000,powerSpecFilt,'EdgeColor','none');
view(2);

xlim([1,sampleNum*5]);
ylim([0,data.range(end)./1000]);

ylim([0 ylimUpper])
xticks(1:sampleNum:sampleNum*5);
grid on

xlabel('Sample number')
ylabel('Range (km)')
title(datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'))

caxis([-80 -25])
colorbar

if saveWaterfall
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'waterfall_',datestr(data.time(startInd),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end

end