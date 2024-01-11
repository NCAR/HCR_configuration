function plotReflectivities(dbzDual,moments,figdir,project,showPlot)

dbzBase=dbzDual(:,:,1);
dbzBaseCheck=dbzDual(:,:,6);
dbzDiff=nan(size(dbzBase));
dbzHigh=nan(size(dbzBase));
dbzLow=nan(size(dbzBase));
dbzHighH=nan(size(dbzBase));
dbzLowH=nan(size(dbzBase));

dbzDiff(:,moments.elevation>0)=-dbzDual(:,moments.elevation>0,7);
dbzDiff(:,moments.elevation<=0)=dbzDual(:,moments.elevation<=0,7);

dbzHigh(:,moments.elevation>0)=dbzDual(:,moments.elevation>0,3);
dbzLow(:,moments.elevation>0)=dbzDual(:,moments.elevation>0,2);
dbzHigh(:,moments.elevation<=0)=dbzDual(:,moments.elevation<=0,2);
dbzLow(:,moments.elevation<=0)=dbzDual(:,moments.elevation<=0,3);

dbzHighH(:,moments.elevation>0)=dbzDual(:,moments.elevation>0,5);
dbzLowH(:,moments.elevation>0)=dbzDual(:,moments.elevation>0,4);
dbzHighH(:,moments.elevation<=0)=dbzDual(:,moments.elevation<=0,4);
dbzLowH(:,moments.elevation<=0)=dbzDual(:,moments.elevation<=0,5);

aslGood=moments.asl(~isnan(dbzHigh))./1000;
ylims=[0,max(aslGood)+0.5];

% lmin=min(dbzLow(:),[],'omitmissing');
% cmin=min([-12,lmin]);
clims=[-30-0.001,20];

f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(4,2,1);

moments.dbz(isnan(moments.dbz))=-99;

colTwo=cat(1,[0,0,0],jet);

hold on
surf(moments.time,moments.asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity time domain (dBZ)')

s2=subplot(4,2,2);

dbzBase(isnan(dbzBase))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzBase,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s2.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity base (dBZ)')

s3=subplot(4,2,3);

dbzHighH(isnan(dbzHighH))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzHighH,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s3.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity vel high (dBZ)')

s4=subplot(4,2,4);

dbzLowH(isnan(dbzLowH))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzLowH,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity vel low (dBZ)')

s5=subplot(4,2,5);

dbzHigh(isnan(dbzHigh))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s5.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity vel high filled (dBZ)')

s6=subplot(4,2,6);

dbzLow(isnan(dbzLow))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s6.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity vel low filled (dBZ)')

s7=subplot(4,2,7);

dbzBaseCheck(isnan(dbzBaseCheck))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzBaseCheck,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s7.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity middle (m s^{-1})')

s8=subplot(4,2,8);

dbzDiff(isnan(dbzDiff))=-99;

hold on
surf(moments.time,moments.asl./1000,dbzDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-15 15]);
s8.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('High minus low (m s^{-1})')

linkaxes([s1 s2 s3 s4 s5 s6 s7 s8],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_dualRefl_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end