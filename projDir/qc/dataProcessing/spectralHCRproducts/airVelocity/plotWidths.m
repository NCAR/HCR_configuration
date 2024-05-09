function plotWidths(moments,momentsSpBasic,momentsSpNoNoise,momentsSpSmooth,momentsSpCorrected,cf,figdir,project,showPlot)
%moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,3];
climsDiff=[-13,13];
colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

%% Figure
f1 = figure('Position',[200 500 2400 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(4,3,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

moments.width(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s1.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width time domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s2=nexttile(2);

momentsSpBasic.width(isnan(momentsSpBasic.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s2.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width spectral domain basic (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpNoNoise.width(isnan(momentsSpNoNoise.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpNoNoise.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width spectral domain noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

cf.WIDTH(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,cf.WIDTH,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width time domain corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);


s5=nexttile(5);

momentsSpCorrected.width(isnan(momentsSpCorrected.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpCorrected.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s5.Colormap=colDiff;
colorbar
grid on
box on
title('Spectrum width spectral domain corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

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

momentsSpBasic.kurt(isnan(momentsSpBasic.kurt))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.kurt,'edgecolor','none');
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