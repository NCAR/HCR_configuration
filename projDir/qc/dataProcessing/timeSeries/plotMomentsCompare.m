function plotMomentsCompare(data,moments,timeBeams,figdir,project,type,ylimUpper,flipYes,showPlot,plotTimes,plotRangeKM)
f1 = figure('Position',[200 500 1800 1300],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet

s1=subplot(4,2,1);

hold on
surf(timeBeams,data.range./1000,moments.powerV,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-110 -40]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Power (dB)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s1.SortMethod='childorder';
end

s2=subplot(4,2,2);

hold on
surf(timeBeams,data.range./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-5 5]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity (m s^{-1})')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s2.SortMethod='childorder';
end

s3=subplot(4,2,3);

hold on
surf(timeBeams,data.range./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-60 20]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Reflectivity (dBZ)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s3.SortMethod='childorder';
end

s4=subplot(4,2,4);

hold on
surf(timeBeams,data.range./1000,moments.width,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([0 4]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Spectrum width (m s^{-1})')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s4.SortMethod='childorder';
end

s5=subplot(4,2,5);

hold on
surf(timeBeams,data.range./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-20 70]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Signal to noise ratio (dB)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s5.SortMethod='childorder';
end

s6=subplot(4,2,6);

hold on
surf(timeBeams,data.range./1000,moments.skew,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-1 1]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Skew (dB)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s6.SortMethod='childorder';
end

s7=subplot(4,2,7);

hold on
surf(timeBeams,data.range./1000,moments.ldr,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-40 10]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('LDR (dB)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s7.SortMethod='childorder';
end

s8=subplot(4,2,8);

hold on
surf(timeBeams,data.range./1000,moments.kurt,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([0 10]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Kurtosis (dB)')

if flipYes
    set(gca, 'YDir','reverse');
end

if ~isempty(plotTimes)
    scatter(plotTimes,plotRangeKM,'ok','LineWidth',1.5);
    s8.SortMethod='childorder';
end

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_moments_',type,'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end