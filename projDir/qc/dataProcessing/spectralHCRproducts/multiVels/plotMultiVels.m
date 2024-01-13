function plotMultiVels(moments,lowShoulderVel,highShoulderVel,figdir,project,showPlot)
%velDual(:,moments.elevation>0,:)=-velDual(:,moments.elevation>0,:);
moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

shoulderLow=nan(size(lowShoulderVel));
shoulderLow(:,moments.elevation<=0)=lowShoulderVel(:,moments.elevation<=0);
shoulderLow(:,moments.elevation>0)=-highShoulderVel(:,moments.elevation>0);
shoulderHigh=nan(size(highShoulderVel));
shoulderHigh(:,moments.elevation<=0)=highShoulderVel(:,moments.elevation<=0);
shoulderHigh(:,moments.elevation>0)=-lowShoulderVel(:,moments.elevation>0);

% velBase=velDual(:,:,1);
% velBaseCheck=velDual(:,:,6);
% velDiff=nan(size(velBase));
% velHigh=nan(size(velBase));
% velLow=nan(size(velBase));
% velHighH=nan(size(velBase));
% velLowH=nan(size(velBase));
% 
% velDiff(:,moments.elevation>0)=-velDual(:,moments.elevation>0,7);
% velDiff(:,moments.elevation<=0)=velDual(:,moments.elevation<=0,7);
% 
% velHigh(:,moments.elevation>0)=velDual(:,moments.elevation>0,3);
% velLow(:,moments.elevation>0)=velDual(:,moments.elevation>0,2);
% velHigh(:,moments.elevation<=0)=velDual(:,moments.elevation<=0,2);
% velLow(:,moments.elevation<=0)=velDual(:,moments.elevation<=0,3);
% 
% velHighH(:,moments.elevation>0)=velDual(:,moments.elevation>0,5);
% velLowH(:,moments.elevation>0)=velDual(:,moments.elevation>0,4);
% velHighH(:,moments.elevation<=0)=velDual(:,moments.elevation<=0,4);
% velLowH(:,moments.elevation<=0)=velDual(:,moments.elevation<=0,5);

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

lmin=min(shoulderLow(:),[],'omitmissing');
cmin=min([-12,lmin]);
clims=[cmin-0.001,abs(cmin)];

f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

s1=subplot(4,2,1);

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

s2=subplot(4,2,2);
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

s3=subplot(4,2,3);

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

s4=subplot(4,2,4);

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

s5=subplot(4,2,5);

% velHigh(isnan(velHigh))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,velHigh,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s5.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity high filled (m s^{-1})')

s6=subplot(4,2,6);

% velLow(isnan(velLow))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,velLow,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s6.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity low filled (m s^{-1})')

s7=subplot(4,2,7);

% velBaseCheck(isnan(velBaseCheck))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,velBaseCheck,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(clims);
% s7.Colormap=colTwo;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Velocity middle (m s^{-1})')

s8=subplot(4,2,8);

% velDiff(isnan(velDiff))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,velDiff,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim([-0.001 8]);
% s8.Colormap=colDiff;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('High minus low (m s^{-1})')

%linkaxes([s1 s2 s3 s4 s5 s6 s7 s8],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiVel_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end