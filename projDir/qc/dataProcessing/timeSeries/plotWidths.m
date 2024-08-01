function plotWidths(moments,momentsSpBasic,momentsSpNoNoise,momentsSpSmooth,momentsSpCorrected,cf,plotTimeAll,figdir,project,showPlot)
%moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsWidth=[0,2];
climsDiff=[-0.8,0.8];
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
title('Time domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

plotRangeInds=[20:20:700];
for kk=1:length(plotTimeAll)
    times=repmat(plotTimeAll(kk),length(plotRangeInds),1);
    alts=momentsSpBasic.asl(plotRangeInds,momentsSpBasic.time==plotTimeAll(kk));
    scatter(times,alts./1000,36,'k','+');
end

s1.SortMethod='childorder';

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
title('Spectral domain basic with noise (m s^{-1})')
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
title('Spectral domain noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

moments.widthCorr(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.widthCorr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Time domain corrected (m s^{-1})')
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
title('Spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s6=nexttile(6);

momentsSpSmooth.width(isnan(momentsSpSmooth.width))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmooth.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s7=nexttile(7);

moments.width(moments.width==-99)=-999;

surf(moments.time,moments.asl./1000,moments.width-moments.widthCorr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s7.Colormap=colTwo;
colorbar
grid on
box on
title('Time domain, raw - corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s8=nexttile(8);
moments.widthCorr(momentsSpCorrected.width==-99)=-99;
moments.widthCorr(moments.widthCorr==-99)=-999;

surf(moments.time,moments.asl./1000,moments.widthCorr-momentsSpCorrected.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s8.Colormap=colTwo;
colorbar
grid on
box on
title('Corrected, time - spectral domain (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s9=nexttile(9);

momentsSpNoNoise.width(momentsSpNoNoise.width==-99)=-999;

surf(moments.time,moments.asl./1000,momentsSpNoNoise.width-momentsSpSmooth.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s9.Colormap=colTwo;
colorbar
grid on
box on
title('Spectral domain, noise removed - filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

% s10=nexttile(10);
% 
% momentsSpNoNoise.width(momentsSpNoNoise.width==-999)=-99;
% 
% surf(moments.time,moments.asl./1000,moments.width-momentsSpNoNoise.width,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(climsDiff);
% s10.Colormap=colTwo;
% colorbar
% grid on
% box on
% title('Time domain - spectral domain noise removed (m s^{-1})')
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);

s10=nexttile(10);

momentsSpNoNoise.velRaw(isnan(momentsSpCorrected.velRaw))=-999;

surf(moments.time,moments.asl./1000,momentsSpNoNoise.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-15,15]);
s10.Colormap=colTwo;
colorbar
grid on
box on
title('Velocity spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s11=nexttile(11);

momentsSpSmooth.width(momentsSpSmooth.width==-99)=-999;

surf(moments.time,moments.asl./1000,momentsSpSmooth.width-momentsSpCorrected.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDiff);
s11.Colormap=colTwo;
colorbar
grid on
box on
title('Spectral domain, filtered - corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s12=nexttile(12);

moments.snr(isnan(moments.vel))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.snr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim([-10,50]);
s12.Colormap=colDiff;
colorbar
grid on
box on
title('SNR (dB)')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_width_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end