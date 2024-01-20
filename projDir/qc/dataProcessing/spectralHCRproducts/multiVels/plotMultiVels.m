function plotMultiVels(moments,shoulderLow,shoulderHigh,velLayers,figdir,project,showPlot)
moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

lmin=min(shoulderLow(:),[],'omitmissing');
cmin=min([-12,lmin]);
clims=[cmin-0.001,abs(cmin)];

f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(3,2,1);

moments.vel(isnan(moments.vel))=-99;

colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

hold on
surf(moments.time,moments.asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity time domain (m s^{-1})')

%s2=subplot(4,2,2);
% 
% velBase(isnan(velBase))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,velBase,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s2.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity base (m s^{-1})')

s3=subplot(3,2,3);

shoulderHigh(isnan(shoulderHigh))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s3.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity high (m s^{-1})')

s4=subplot(3,2,4);

shoulderLow(isnan(shoulderLow))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity low (m s^{-1})')

s5=subplot(3,2,5);

highLayer=velLayers(:,:,2);
highLayer(isnan(highLayer))=-99;

hold on
surf(moments.time,moments.asl./1000,highLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s5.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles high (m s^{-1})')

s6=subplot(3,2,6);

lowLayer=velLayers(:,:,1);
lowLayer(isnan(lowLayer))=-99;

hold on
surf(moments.time,moments.asl./1000,lowLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s6.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles low (m s^{-1})')

linkaxes([s1 s3 s4 s5 s6],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiVel_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end