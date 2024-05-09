function plotAirMotion(moments,momentsSp,cf,shoulderLow,shoulderHigh,peakLow,peakHigh,aircraft,figdir,project,showPlot,plotTimeAll)
%moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=momentsSp.asl(~isnan(momentsSp.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

clims=[-13,13];
colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

dualPartDiff=peakHigh-peakLow;

%% Set up echo type and melting layer
colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    0.7,0,0];

stratConvPlot=cf.ECHO_TYPE_2D;
stratConvPlot(stratConvPlot==14)=1;
stratConvPlot(stratConvPlot==16)=2;
stratConvPlot(stratConvPlot==18)=3;
stratConvPlot(stratConvPlot==25)=4;
stratConvPlot(stratConvPlot==30)=5;
stratConvPlot(stratConvPlot==32)=6;
stratConvPlot(stratConvPlot==34)=7;
stratConvPlot(stratConvPlot==36)=8;
stratConvPlot(stratConvPlot==38)=9;

% Melting layer
meltPlot=cf.MELTING_LAYER;
meltPlot(meltPlot==9)=3;
meltPlot(meltPlot==11)=2;
meltPlot(meltPlot==19)=1;
meltPlot(meltPlot==21)=0;

%% Figure
f1 = figure('Position',[200 500 2400 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

t = tiledlayout(4,3,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

momentsSp.velRaw(isnan(momentsSp.velRaw))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s1.Colormap=colTwo;
ylim(ylims);
xlim([momentsSp.time(1),momentsSp.time(end)]);
colorbar
grid on
box on
title('Velocity (m s^{-1})')

plotRangeInds=[20:20:700];
for kk=1:length(plotTimeAll)
    times=repmat(plotTimeAll(kk),length(plotRangeInds),1);
    alts=momentsSp.asl(plotRangeInds,momentsSp.time==plotTimeAll(kk));
    scatter(times,alts./1000,36,'b','x');
end

scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
set(gca,'clim',clims);

s1.SortMethod='childorder';

s2=nexttile(2);

momentsSp.width(isnan(momentsSp.width))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,3]);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width (m s^{-1})')

ax2 = axes(t);
ax2.Layout.Tile=2;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax2.Colormap=colTwo;
ax2.Visible = 'off';
ax2.CLim=clims;
linkaxes([s2,ax2]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s2.SortMethod='childorder';

s3=nexttile(3);

moments.dbz(isnan(moments.vel))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.dbz,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-40,30]);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Reflectivity (dBZ)')

ax3=axes(t);
ax3.Layout.Tile=3;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax3.Colormap=colTwo;
ax3.Visible = 'off';
ax3.CLim=clims;
linkaxes([s3,ax3]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3.SortMethod='childorder';

s4=nexttile(4);

shoulderHigh(isnan(shoulderHigh))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderHigh,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s4.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity high (m s^{-1})')

scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
set(gca,'clim',clim);

s4.SortMethod='childorder';

s5=nexttile(5);

momentsSp.skew(isnan(momentsSp.skew))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-3,3]);
s5.Colormap=colTwo;
colorbar
grid on
box on
title('Skew (dB)')

ax5=axes(t);
ax5.Layout.Tile=5;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax5.Colormap=colTwo;
ax5.Visible = 'off';
ax5.CLim=clims;
linkaxes([s5,ax5]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s5.SortMethod='childorder';

s6=nexttile(6);

moments.snr(isnan(moments.vel))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-20,60]);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('SNR (dB)')

ax6=axes(t);
ax6.Layout.Tile=6;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax6.Colormap=colTwo;
ax6.Visible = 'off';
ax6.CLim=clims;
linkaxes([s6,ax6]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s6.SortMethod='childorder';

s7=nexttile(7);

shoulderLow(isnan(shoulderLow))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s7.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity low (m s^{-1})')

scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
set(gca,'clim',clim);

s8=nexttile(8);

momentsSp.kurt(isnan(momentsSp.kurt))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-6,6]);
s8.Colormap=colTwo;
colorbar
grid on
box on
title('Kurtosis (dB)')

ax8=axes(t);
ax8.Layout.Tile=8;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax8.Colormap=colTwo;
ax8.Visible = 'off';
ax8.CLim=clims;
linkaxes([s8,ax8]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s8.SortMethod='childorder';

s9=nexttile(9);

cf.CONVECTIVITY(isnan(cf.ECHO_TYPE_2D))=-99;

hold on
surf(cf.time,cf.asl./1000,cf.CONVECTIVITY,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,1]);
s9.Colormap=colDiff;
colorbar
grid on
box on
title('Convectivity')

ax9=axes(t);
ax9.Layout.Tile=9;
scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
ax9.Colormap=colTwo;
ax9.Visible = 'off';
ax9.CLim=clims;
linkaxes([s9,ax9]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s9.SortMethod='childorder';

s10=nexttile(10);

peakLow(isnan(peakLow))=-99;

hold on
surf(moments.time,moments.asl./1000,peakLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s10.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles low (m s^{-1})')

scatter(aircraft.Time,aircraft.Alt./1000,20,-aircraft.Vel,'filled');
set(gca,'clim',clim);

s10.SortMethod='childorder';

% s10=nexttile(10);
% 
% dualPartDiff(isnan(dualPartDiff))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,dualPartDiff,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim([0,6]);
% s10.Colormap=colDiff;
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);
% colorbar
% grid on
% box on
% title('Dual particles high-low (m s^{-1})')

s11=nexttile(11);

hold on;
surf(cf.time,cf.asl./1000,meltPlot,'edgecolor','none');
view(2);
s11.Colormap=[0,1,1;0.5,0.5,0.5;0,0,0;1,0,1];
clim([-0.5,3.5]);
cb=colorbar;
cb.Ticks=[0,1,2,3];
cb.TickLabels={'Cold','Melting cold','Melting warm','Warm'};
ylim(ylims);
ylabel('Altitude (km)');
xlim([moments.time(1),moments.time(end)]);
title('Melting layer')
grid on
box on

s12=nexttile(12);

hold on
surf(cf.time,cf.asl./1000,stratConvPlot,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0.5 9.5]);
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
s12.Colormap=colmapSC;
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels={'StratLow','StratMid','StratHigh','Mixed',...
    'Conv','ConvElev','ConvShallow','ConvMid','ConvDeep'};
grid on
box on
title('Stratiform/convective echo type')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_airMotion_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end