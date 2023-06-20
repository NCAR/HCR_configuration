% Analyze HCR clouds

clear all;
close all;

project='cset'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

% figdir=['/scr/snow1/rsfdata/projects/otrec/hcr/qc2/cfradial/final2/10hz/plots/testHourly/'];
figdir=['/home/romatsch/plots/HCR/meltingLayer/hourly/',project,'/',freqData,'/dropSondeComp/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

ylimits=[-0.2 7];
        
compAlts=readtable([figdir,project,'_meltLayer_dropsonde.txt']);

%% Plot comparison scatter plot
close all

fig2=figure('DefaultAxesFontSize',11,'position',[100,1300,700,650]);
hold on

maxLim=ceil(max(max(table2array(compAlts(:,2:end))))/1000);
minLim=floor(min(min(table2array(compAlts(:,2:end))))/1000);
offset=(maxLim-minLim)*0.04;
l0=plot([minLim maxLim],[minLim maxLim],'-r');
xticks([minLim:0.5:maxLim]);
yticks([minLim:0.5:maxLim]);

% Regression lines
fitOrth1=gmregress(compAlts.zeroDegAlt./1000,compAlts.sondeAlt./1000,1);
fitAll1=[fitOrth1(2) fitOrth1(1)];
xFitD =minLim:0.1:maxLim;
yFitD1 = polyval(fitAll1, xFitD);
l1=plot(xFitD, yFitD1,'-k','linewidth',1.5);
fitOrth2=gmregress(compAlts.meltAltEst./1000,compAlts.sondeAlt./1000,1);
fitAll2=[fitOrth2(2) fitOrth2(1)];
yFitD2 = polyval(fitAll2, xFitD);
l2=plot(xFitD, yFitD2,'-g','linewidth',1.5);
fitOrth3=gmregress(compAlts.meltAltInt./1000,compAlts.sondeAlt./1000,1);
fitAll3=[fitOrth3(2) fitOrth3(1)];
yFitD3 = polyval(fitAll3, xFitD);
l3=plot(xFitD, yFitD3,'-c','linewidth',1.5);
fitOrth4=gmregress(compAlts.meltAltMeas./1000,compAlts.sondeAlt./1000,1);
fitAll4=[fitOrth4(2) fitOrth4(1)];
yFitD4 = polyval(fitAll4, xFitD);
l4=plot(xFitD, yFitD4,'-b','linewidth',1.5);

scatter(compAlts.zeroDegAlt./1000,compAlts.sondeAlt./1000,30,'k','filled');
scatter(compAlts.meltAltEst./1000,compAlts.sondeAlt./1000,30,'g','filled');
scatter(compAlts.meltAltInt./1000,compAlts.sondeAlt./1000,30,'c','filled');
scatter(compAlts.meltAltMeas./1000,compAlts.sondeAlt./1000,30,'b','filled');
xlabel('HCR altitude (km)');
ylabel('Dropsonde altitude (km)');

title(['Melting layer altitude'])
grid on
axis equal

xlim([minLim maxLim]);
ylim([minLim maxLim]);

measDiff=mean((compAlts.meltAltMeas-compAlts.sondeAlt)./1000,'omitnan');
measStd=std((compAlts.meltAltMeas-compAlts.sondeAlt)./1000,'omitnan');
intDiff=mean((compAlts.meltAltInt-compAlts.sondeAlt)./1000,'omitnan');
intStd=std((compAlts.meltAltInt-compAlts.sondeAlt)./1000,'omitnan');
estDiff=mean((compAlts.meltAltEst-compAlts.sondeAlt)./1000,'omitnan');
estStd=std((compAlts.meltAltEst-compAlts.sondeAlt)./1000,'omitnan');
zeroDiff=mean((compAlts.zeroDegAlt-compAlts.sondeAlt)./1000,'omitnan');
zeroStd=std((compAlts.zeroDegAlt-compAlts.sondeAlt)./1000,'omitnan');

text(minLim+offset,maxLim-offset,['Measured - dropsonde: ',num2str(measDiff,2), ...
    ' +/- ',num2str(measStd,2),' km'],'fontsize',12);
text(minLim+offset,maxLim-2*offset,['Interpolated - dropsonde: ',num2str(intDiff,2), ...
    ' +/- ',num2str(intStd,2),' km'],'fontsize',12);
text(minLim+offset,maxLim-3*offset,['Estimated - dropsonde: ',num2str(estDiff,2), ...
    ' +/- ',num2str(estStd,2),' km'],'fontsize',12);
text(minLim+offset,maxLim-4*offset,['Zero deg - dropsonde: ',num2str(zeroDiff,2), ...
    ' +/- ',num2str(zeroStd,2),' km'],'fontsize',12);

leg1=legend([l4,l3,l2,l1],{'Measured','Interpolated','Estimated','Zero deg'},'Location','southeast');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_meltLayer_dropsonde.png'],'-dpng','-r0');
