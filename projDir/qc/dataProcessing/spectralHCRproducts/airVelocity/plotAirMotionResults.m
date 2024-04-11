function plotAirMotionResults(moments,airMotion,TTcompAC,TTcloud,TTprecip,figdir,project,showPlot)
%moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

clims=[-13,13];
colTwo=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1200 650],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
s1=nexttile([1,2]);

airMotion(isnan(airMotion))=-99;

hold on
surf(moments.time,moments.asl./1000,airMotion,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Air motion (m s^{-1})')

scatter(TTcompAC.Time,TTcompAC.aircraftAlt./1000,20,-TTcompAC.aircraftVel,'filled');
set(gca,'clim',clims);

s1.SortMethod='childorder';

s2=nexttile(4);

scatter(-TTcloud.aircraftVel,TTcloud.totVel,'filled');
xlimsIn=s2.XLim;
ylimsIn=s2.YLim;
lims=[min([xlimsIn(1),ylimsIn(1)]),max([xlimsIn(2),ylimsIn(2)])];
xlim(lims);
ylim(lims);
title('Cloud')
xlabel('Velocity from aircraft (m s^{-1})');
ylabel('Velocity from radar (m s^{-1})');
grid on
box on

s3=nexttile(5);

scatter(-TTprecip.aircraftVel,TTprecip.lowVel,'filled');
xlimsIn=s3.XLim;
ylimsIn=s3.YLim;
lims=[min([xlimsIn(1),ylimsIn(1)]),max([xlimsIn(2),ylimsIn(2)])];
xlim(lims);
ylim(lims);
title('Precip')
xlabel('Velocity from aircraft (m s^{-1})');
ylabel('Velocity from radar (m s^{-1})');
grid on
box on

s4=nexttile(6);

scatter(-TTcompAC.aircraftVel,TTcompAC.velCombined,'filled');
xlimsIn=s4.XLim;
ylimsIn=s4.YLim;
lims=[min([xlimsIn(1),ylimsIn(1)]),max([xlimsIn(2),ylimsIn(2)])];
xlim(lims);
ylim(lims);
title('Combined')
xlabel('Velocity from aircraft (m s^{-1})');
ylabel('Velocity from radar (m s^{-1})');
grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_airMotionResults_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end