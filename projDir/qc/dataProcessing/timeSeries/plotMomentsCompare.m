function plotMomentsCompare(moments,figdir,project,type,ylimUpper,showPlot,plotTimes,plotRangeKM)
f1 = figure('Position',[200 500 1800 1300],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

s1=subplot(4,2,1);

hold on
surf(moments.time,moments.asl./1000,moments.powerV,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-110 -40]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Power (dB)')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s1.SortMethod='childorder';
end

s2=subplot(4,2,2);

hold on
surf(moments.time,moments.asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-5 5]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity (m s^{-1})')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s2.SortMethod='childorder';
end

s3=subplot(4,2,3);

hold on
surf(moments.time,moments.asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-60 20]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity (dBZ)')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s3.SortMethod='childorder';
end

s4=subplot(4,2,4);

hold on
surf(moments.time,moments.asl./1000,moments.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0 4]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Spectrum width (m s^{-1})')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s4.SortMethod='childorder';
end

s5=subplot(4,2,5);

hold on
surf(moments.time,moments.asl./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-20 70]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Signal to noise ratio (dB)')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s5.SortMethod='childorder';
end

s6=subplot(4,2,6);

hold on
surf(moments.time,moments.asl./1000,moments.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-1 1]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Skew (dB)')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s6.SortMethod='childorder';
end

if ~isempty(find(~isnan(moments.ldr)))
    s7=subplot(4,2,7);

    hold on
    surf(moments.time,moments.asl./1000,moments.ldr,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim([-40 10]);
    ylim([0 ylimUpper]);
    xlim([moments.time(1),moments.time(end)]);
    colorbar
    grid on
    box on
    title('LDR (dB)')

    if ~isempty(plotTimes)
        scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
        s7.SortMethod='childorder';
    end
end

s8=subplot(4,2,8);

hold on
surf(moments.time,moments.asl./1000,moments.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0 10]);
ylim([0 ylimUpper]);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Kurtosis (dB)')

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s8.SortMethod='childorder';
end

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_moments_',type,'_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end