function plotAllMoments(cf,momentsSpec,figdir,project,showPlot)

momentsSpec.velRaw(:,momentsSpec.elevation>0)=-momentsSpec.velRaw(:,momentsSpec.elevation>0);
momentsSpec.skew(:,momentsSpec.elevation>0)=-momentsSpec.skew(:,momentsSpec.elevation>0);

aslGood=momentsSpec.asl(~isnan(momentsSpec.velRaw))./1000;
ylims=[0,max(aslGood)+0.5];

climsDbz=[-40,30];
climsLdr=[-35,-5];
climsVel=[-15,15];
climsWidth=[0,3];
climsSkew=[-3,3];
climsKurt=[-6,6];

col1=cat(1,[0,0,0],jet);
col2=cat(1,[0,0,0],velCols);

%% Figure
f1 = figure('Position',[200 500 1400 1000],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

cf.DBZ(isnan(cf.VEL_MASKED))=-999;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,cf.DBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsDbz);
s1.Colormap=col1;
colorbar
grid on
box on
title('Time domain reflectivity (dBZ)')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

if isfield(cf,'LDR')
    s2=nexttile(2);

    cf.LDR(isnan(cf.VEL_MASKED))=-999;
    cf.LDR(isnan(cf.LDR))=-999;

    hold on
    surf(momentsSpec.time,momentsSpec.asl./1000,cf.LDR,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    clim(climsLdr);
    s2.Colormap=col1;
    colorbar
    grid on
    box on
    title('Time domain linear depolarization ratio (dB)')
    ylim(ylims);
    xlim([momentsSpec.time(1),momentsSpec.time(end)]);
end

s3=nexttile(3);

momentsSpec.velRaw(isnan(momentsSpec.velRaw))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.velRaw,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsVel);
s3.Colormap=col2;
colorbar
grid on
box on
title('Spectral domain velocity (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s4=nexttile(4);

momentsSpec.width(isnan(momentsSpec.width))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.width,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsWidth);
s4.Colormap=col1;
colorbar
grid on
box on
title('Spectral domain width (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s5=nexttile(5);

momentsSpec.skew(isnan(momentsSpec.skew))=-99;
momentsSpec.skew(isinf(momentsSpec.skew))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.skew,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s5.Colormap=col2;
colorbar
grid on
box on
title('Spectral domain skewness (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);

s6=nexttile(6);

momentsSpec.kurt(isnan(momentsSpec.kurt))=-99;
momentsSpec.kurt(isinf(momentsSpec.kurt))=-99;

hold on
surf(momentsSpec.time,momentsSpec.asl./1000,momentsSpec.kurt,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s6.Colormap=col2;
colorbar
grid on
box on
title('Spectral domain kurtosis (m s^{-1})')
ylim(ylims);
xlim([momentsSpec.time(1),momentsSpec.time(end)]);


set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_moments_',datestr(momentsSpec.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(momentsSpec.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end