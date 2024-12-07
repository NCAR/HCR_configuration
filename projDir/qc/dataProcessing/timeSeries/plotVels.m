function plotVels(moments,momentsSpBasic,momentsSpBasicRMnoise,momentsSpSmooth,momentsSpSmoothCorr,cf,plotTimeAll,figdir,project,showPlot)
moments.vel(:,moments.elevation>0,:)=-moments.vel(:,moments.elevation>0,:);
momentsSpBasic.velRaw(:,moments.elevation>0,:)=-momentsSpBasic.velRaw(:,moments.elevation>0,:);
momentsSpBasicRMnoise.velRaw(:,moments.elevation>0,:)=-momentsSpBasicRMnoise.velRaw(:,moments.elevation>0,:);
momentsSpSmooth.velRaw(:,moments.elevation>0,:)=-momentsSpSmooth.velRaw(:,moments.elevation>0,:);
momentsSpSmoothCorr.velRaw(:,moments.elevation>0,:)=-momentsSpSmoothCorr.velRaw(:,moments.elevation>0,:);

aslGood=momentsSpBasic.asl(~isnan(momentsSpBasic.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsVel=[-15,15];
colDiff=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1600 1250],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

moments.vel(isnan(cf.VEL_MASKED))=-99;

hold on
surf(moments.time,moments.asl./1000,moments.vel,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
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

% s2=nexttile(2);
% 
% moments.widthCorr(isnan(cf.VEL_MASKED))=-99;
% 
% hold on
% surf(moments.time,moments.asl./1000,moments.widthCorr,'edgecolor','none');
% view(2);
% ylabel('Altitude (km)');
% clim(climsVel);
% s2.Colormap=colDiff;
% colorbar
% grid on
% box on
% title('Time domain corrected (m s^{-1})')
% ylim(ylims);
% xlim([moments.time(1),moments.time(end)]);

s3=nexttile(3);

momentsSpBasic.velRaw(isnan(momentsSpBasic.velRaw))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasic.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s3.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s4=nexttile(4);

momentsSpBasicRMnoise.velRaw(isnan(momentsSpBasicRMnoise.velRaw))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpBasicRMnoise.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s4.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain raw noise removed (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s5=nexttile(5);

momentsSpSmooth.velRaw(isnan(momentsSpSmooth.velRaw))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmooth.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s5.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

s6=nexttile(6);

momentsSpSmoothCorr.velRaw(isnan(momentsSpSmoothCorr.velRaw))=-99;

hold on
surf(momentsSpBasic.time,momentsSpBasic.asl./1000,momentsSpSmoothCorr.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('Spectral domain filtered and corrected (m s^{-1})')
ylim(ylims);
xlim([moments.time(1),moments.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_vel_',datestr(moments.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(moments.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end