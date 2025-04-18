clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

fo=load([figdir,'forReviewers/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);
fm1=load([figdir,'forReviewers/filterMinus1/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);
fm2=load([figdir,'forReviewers/filterMinus2/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);
fp1=load([figdir,'forReviewers/filterPlus1/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);
fp2=load([figdir,'forReviewers/filterPlus2/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);
filC=load([figdir,'forReviewers/specData/noreaster/noreaster_filterCorr.mat']);
snrIn=load([figdir,'widthNoreaster.mat']);

%% Match data
startTime=datetime(2015,2,2,12,55,0);
endTime=snrIn.momentsTime.time(end);

timeInds=find(fo.momentsSpecParams.time>=startTime & fo.momentsSpecParams.time<=endTime);
timeInds2=find(snrIn.momentsTime.time>=startTime & snrIn.momentsTime.time<=endTime);

climsSkew=[-3,3];
climsKurt=[-6,6];

col1=jet;
col2=velCols;

xlims=[datetime(2015,2,2,12,55,26),datetime(2015,2,2,12,56,59)];
ylims=[1.5,10];
%% Figure
close all

f1 = figure('Position',[200 500 1100 1200],'DefaultAxesFontSize',12,'visible','on');

t = tiledlayout(6,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

surf(fm2.momentsSpecParams.time(timeInds),fm2.momentsSpecParams.asl(:,timeInds)./1000,fm2.momentsSpecParams.skew(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s1.Colormap=col2;
colorbar
grid on
box on
title('(a) Optimum TRV-2; Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s2=nexttile(2);

surf(fm2.momentsSpecParams.time(timeInds),fm2.momentsSpecParams.asl(:,timeInds)./1000,fm2.momentsSpecParams.kurt(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s2.Colormap=col2;
colorbar
grid on
box on
title('(b) Optimum TRV-2; Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);

s3=nexttile(3);

surf(fm1.momentsSpecParams.time(timeInds),fm1.momentsSpecParams.asl(:,timeInds)./1000,fm1.momentsSpecParams.skew(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s3.Colormap=col2;
colorbar
grid on
box on
title('(c) Optimum TRV-1; Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s4=nexttile(4);

surf(fm1.momentsSpecParams.time(timeInds),fm1.momentsSpecParams.asl(:,timeInds)./1000,fm1.momentsSpecParams.kurt(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s4.Colormap=col2;
colorbar
grid on
box on
title('(d) Optimum TRV-1; Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);

s5=nexttile(5);

surf(fo.momentsSpecParams.time(timeInds),fo.momentsSpecParams.asl(:,timeInds)./1000,fo.momentsSpecParams.skew(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s5.Colormap=col2;
colorbar
grid on
box on
title('(e) Optimum TRV; Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s6=nexttile(6);

surf(fo.momentsSpecParams.time(timeInds),fo.momentsSpecParams.asl(:,timeInds)./1000,fo.momentsSpecParams.kurt(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s6.Colormap=col2;
colorbar
grid on
box on
title('(f) Optimum TRV; Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);

s7=nexttile(7);

surf(fp1.momentsSpecParams.time(timeInds),fp1.momentsSpecParams.asl(:,timeInds)./1000,fp1.momentsSpecParams.skew(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s7.Colormap=col2;
colorbar
grid on
box on
title('(g) Optimum TRV+1; Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s8=nexttile(8);

surf(fp1.momentsSpecParams.time(timeInds),fp1.momentsSpecParams.asl(:,timeInds)./1000,fp1.momentsSpecParams.kurt(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s8.Colormap=col2;
colorbar
grid on
box on
title('(h) Optimum TRV+1; Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);

s9=nexttile(9);

surf(fp2.momentsSpecParams.time(timeInds),fp2.momentsSpecParams.asl(:,timeInds)./1000,fp2.momentsSpecParams.skew(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsSkew);
s9.Colormap=col2;
colorbar
grid on
box on
title('(i) Optimum TRV+2; Skewness (m s^{-1})')
ylim(ylims);
xlim(xlims);

s10=nexttile(10);

surf(fp2.momentsSpecParams.time(timeInds),fp2.momentsSpecParams.asl(:,timeInds)./1000,fp2.momentsSpecParams.kurt(:,timeInds),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
clim(climsKurt);
s10.Colormap=col2;
colorbar
grid on
box on
title('(j) Optimum TRV+2; Kurtosis (m s^{-1})')
ylim(ylims);
xlim(xlims);

s11=nexttile(11);
hold on
scatter(fp2.momentsSpecParams.time(timeInds),filC.corrFactor(timeInds),'filled');
ylim([0.55,0.65]);
ylabel('BBS (m s^{-1})')

yyaxis right
scatter(fp2.momentsSpecParams.time(timeInds),filC.velTestWind(timeInds),'filled');
box on
xlim(xlims);
ylabel('Aircraft speed (m s^{-1})')

title('(k) BBS and aircraft speed')

s12=nexttile(12);
snr=snrIn.momentsTime.snr(:,timeInds2);
snr(isnan(fp2.momentsSpecParams.kurt(:,timeInds)))=nan;
surf(fp2.momentsSpecParams.time(timeInds),fp2.momentsSpecParams.asl(:,timeInds)./1000,snr,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
%clim(climsKurt);
s12.Colormap=col1;
colorbar
grid on
box on
title('(l) Signal to noise ratio (dB)')
ylim(ylims);
xlim(xlims);

set(gcf,'PaperPositionMode','auto')
%print(f1,[figdir,'filterExample.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'filterExample.png'],'Resolution','300');