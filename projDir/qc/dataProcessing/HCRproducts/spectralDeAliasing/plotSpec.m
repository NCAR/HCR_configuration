function plotSpec(data,sampleNum,startInd,powerSpecLarge,ylimUpper,powerSpecFilt,powerSpecMed,powerSpecMed2,plotGates,showPlot,figdir)
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

caxis([-100 25])
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

caxis([-100 25])
colorbar

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'waterfall_',datestr(data.time(startInd),'yyyymmdd_HHMMSS')],'-dpng','-r0');

%% Plot spectra

if plotGates
    f1 = figure('Position',[100 500 600 400],'DefaultAxesFontSize',12,'visible',showPlot);

    for kk=1:770

        scatter(1:sampleNum*5,powerSpecLarge(kk,:),'filled');
        hold on
        scatter(1:sampleNum*5,powerSpecFilt(kk,:),'filled');
        plot(1:sampleNum*5,powerSpecMed(kk,:),'-k','linewidth',1.5);
        plot(1:sampleNum*5,powerSpecMed2(kk,:),'-c','linewidth',1.5);
        xlabel('Sample number')
        ylabel('Power (dB)')

        ylim([-100 0]);
        hold off

        title([datestr(data.time(startInd),'yyyy-mm-dd HH:MM:SS'),' range ',num2str(data.range(kk)./1000,2),' km'])

        %print(f1,[figdir,'spectrum_',datestr(data.time(startInd),'yyyymmdd_HHMMSS'),'_range_',num2str(data.range(plotRangeInd)./1000),'km'],'-dpng','-r0');

    end
end
end