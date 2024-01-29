function plotMultiRefs(moments,shoulderLowDbz,shoulderHighDbz,dbzLayers,figdir,project,showPlot)

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

% lmin=min(shoulderLowPow(:),[],'omitmissing');
% cmin=min([-12,lmin]);
% clims=[cmin-0.001,abs(cmin)];
clims=[-20,25];

f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(3,2,1);

moments.dbz(isnan(moments.vel))=-99;

%colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

hold on
surf(moments.time,moments.asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity time domain (dBZ)')

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
% title('Reflectivity base (dBZ)')

s3=subplot(3,2,3);

shoulderHighDbz(isnan(shoulderHighDbz))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderHighDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s3.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity high (dBZ)')

s4=subplot(3,2,4);

shoulderLowDbz(isnan(shoulderLowDbz))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderLowDbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Reflectivity low (dBZ)')

s5=subplot(3,2,5);

highLayer=dbzLayers(:,:,2);
highLayer(isnan(highLayer))=-99;

hold on
surf(moments.time,moments.asl./1000,highLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s5.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles high (dBZ)')

s6=subplot(3,2,6);

lowLayer=dbzLayers(:,:,1);
lowLayer(isnan(lowLayer))=-99;

hold on
surf(moments.time,moments.asl./1000,lowLayer,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s6.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles low (dBZ)')

linkaxes([s1 s3 s4 s5 s6],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiRef_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end