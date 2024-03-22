function plotMultiVels(moments,momentsSp,shoulderLow,shoulderHigh,peakLow,peakHigh,peakPowLow,peakPowHigh,figdir,project,showPlot,plotTimeAll)
moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=moments.asl(~isnan(moments.vel))./1000;
ylims=[0,max(aslGood)+0.5];

% lmin=min(shoulderLow(:),[],'omitmissing');
% cmin=min([-12,lmin]);
% clims=[cmin-0.001,abs(cmin)];
clims=[-13,13];

f1 = figure('Position',[200 500 2400 1250],'DefaultAxesFontSize',12,'visible',showPlot);

colormap(velCols);

t = tiledlayout(4,3,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

moments.vel(isnan(moments.vel))=-99;

colTwo=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet);

dualPartDiff=peakHigh-peakLow;
dualPartPowDiff=peakPowHigh-peakPowLow;
shoulderDiff=shoulderHigh-shoulderLow;


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

plotRangeInds=[20:20:700];
for kk=1:length(plotTimeAll)
    times=repmat(plotTimeAll(kk),length(plotRangeInds),1);
    alts=moments.asl(plotRangeInds,moments.time==plotTimeAll(kk));
scatter(times,alts./1000,36,'b','x');
end
s1.SortMethod='childorder';

s2=nexttile(2);

momentsSp.velRaw(isnan(momentsSp.velRaw))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s2.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity spectral domain (m s^{-1})')

s3=nexttile(3);

momentsSp.width(isnan(momentsSp.width))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,3]);
s3.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Spectrum width (m s^{-1})')

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

s5=nexttile(5);

peakHigh(isnan(peakHigh))=-99;

hold on
surf(moments.time,moments.asl./1000,peakHigh,'edgecolor','none');
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

s6=nexttile(6);

momentsSp.skew(isnan(momentsSp.skew))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-2,2]);
s6.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Skew (dB)')

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

s8=nexttile(8);

peakLow(isnan(peakLow))=-99;

hold on
surf(moments.time,moments.asl./1000,peakLow,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(clims);
s8.Colormap=colTwo;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles low (m s^{-1})')

s9=nexttile(9);

momentsSp.kurt(isnan(momentsSp.kurt))=-99;

hold on
surf(momentsSp.time,momentsSp.asl./1000,momentsSp.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,8]);
s9.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Kurtosis (dB)')

s10=nexttile(10);

shoulderDiff(isnan(shoulderDiff))=-99;

hold on
surf(moments.time,moments.asl./1000,shoulderDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,12]);
s10.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Velocity high-low (m s^{-1})')

s11=nexttile(11);

dualPartDiff(isnan(dualPartDiff))=-99;

hold on
surf(moments.time,moments.asl./1000,dualPartDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([0,6]);
s11.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles high-low (m s^{-1})')

s12=nexttile(12);

dualPartPowDiff(isnan(dualPartPowDiff))=-99;

hold on
surf(moments.time,moments.asl./1000,dualPartPowDiff,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-15,15]);
s12.Colormap=colDiff;
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);
colorbar
grid on
box on
title('Dual particles power high-low (dB)')

linkaxes([s1 s4 s5 s6 s7 s8 s9 s10 s11 s12],'xy')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_multiVel_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end