function plotVelocities(data,velDual,moments,timeBeams,figdir,project,ylimUpper,flipYes,showPlot)
velMin=min(velDual,[],3,'omitmissing');
velMax=max(velDual,[],3,'omitmissing');

velMin(isnan(velMin) & ~isnan(velMax))=velMax(isnan(velMin) & ~isnan(velMax));
velMax(isnan(velMax) & ~isnan(velMin))=velMin(isnan(velMax) & ~isnan(velMin));

f1 = figure('Position',[200 500 1800 900],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(2,2,1);

hold on
surf(timeBeams,data.range./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-8 8]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity time domain (m s^{-1})')

if flipYes
    set(gca, 'YDir','reverse');
end

s2=subplot(2,2,2);

hold on
surf(timeBeams,data.range./1000,velMin,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-8 8]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity minimum (m s^{-1})')

if flipYes
    set(gca, 'YDir','reverse');
end

s3=subplot(2,2,3);

overlap=sum(~isnan(velDual),3);
plotOverlap=double(overlap>1);
velOverlap=moments.vel;
velOverlap(plotOverlap==0)=nan;
velOverlap(~isnan(moments.vel) & isnan(velOverlap))=-99;

colTwo=cat(1,[0,1,0],velCols);

hold on
surf(timeBeams,data.range./1000,velOverlap,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-8.001 8]);
s3.Colormap=colTwo;
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Regions with two particles')

if flipYes
    set(gca, 'YDir','reverse');
end

s4=subplot(2,2,4);

hold on
surf(timeBeams,data.range./1000,velMax,'edgecolor','none');
view(2);
ylabel('Range (km)');
caxis([-8 8]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title('Velocity maximum (m s^{-1})')

if flipYes
    set(gca, 'YDir','reverse');
end

linkaxes([s1 s2 s3 s4],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_vels_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end