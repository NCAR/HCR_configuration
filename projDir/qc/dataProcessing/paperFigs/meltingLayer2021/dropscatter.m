% Analyze HCR clouds

clear all;
close all;

freqData='10hz'; % 10hz, 100hz, 2hz, or combined

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));
figdir='/home/romatsch/plots/HCR/meltingLayer/paper/';

%% Read data
project='cset'; %socrates, aristo, cset
indir=['/home/romatsch/plots/HCR/meltingLayer/dropsonde/',project,'/'];
        
compAltsC=readtable([indir,project,'_meltLayer_dropsonde.txt']);

project='socrates'; %socrates, aristo, cset
indir=['/home/romatsch/plots/HCR/meltingLayer/dropsonde/',project,'/'];
        
compAltsS=readtable([indir,project,'_meltLayer_dropsonde.txt']);

project='otrec'; %socrates, aristo, cset
indir=['/home/romatsch/plots/HCR/meltingLayer/dropsonde/',project,'/'];
        
compAltsO=readtable([indir,project,'_meltLayer_dropsonde.txt']);

%% Plot comparison scatter plot
close all
wi=5;
hi=8;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

formatSpec = '%.2f';

set(fig1,'color','w');

ax1=subplot(3,1,1);
hold on;
maxLim=ceil(max(max(table2array(compAltsC(:,2:end))))/1000);
minLim=floor(min(min(table2array(compAltsC(:,2:end))))/1000);
offset=(maxLim-minLim)*0.08;
l0=plot([minLim maxLim],[minLim maxLim],'-r');
xticks([minLim:0.5:maxLim]);
yticks([minLim:0.5:maxLim]);

% Regression lines
fitOrth1=gmregress(compAltsC.zeroDegAlt./1000,compAltsC.sondeAlt./1000,1);
fitAll1=[fitOrth1(2) fitOrth1(1)];
xFitD =minLim:0.1:maxLim;
yFitD1 = polyval(fitAll1, xFitD);
l1=plot(xFitD, yFitD1,'-k','linewidth',1.5);
fitOrth2=gmregress(compAltsC.meltAltEst./1000,compAltsC.sondeAlt./1000,1);
fitAll2=[fitOrth2(2) fitOrth2(1)];
yFitD2 = polyval(fitAll2, xFitD);
l2=plot(xFitD, yFitD2,'color',[0.2 0.6 0.04],'linewidth',1.5);
fitOrth3=gmregress(compAltsC.meltAltInt./1000,compAltsC.sondeAlt./1000,1);
fitAll3=[fitOrth3(2) fitOrth3(1)];
yFitD3 = polyval(fitAll3, xFitD);
l3=plot(xFitD, yFitD3,'-c','linewidth',1.5);
fitOrth4=gmregress(compAltsC.meltAltMeas./1000,compAltsC.sondeAlt./1000,1);
fitAll4=[fitOrth4(2) fitOrth4(1)];
yFitD4 = polyval(fitAll4, xFitD);
l4=plot(xFitD, yFitD4,'-b','linewidth',1.5);

scatter(compAltsC.zeroDegAlt./1000,compAltsC.sondeAlt./1000,30,'k','filled','MarkerEdgeColor','k');
scatter(compAltsC.meltAltEst./1000,compAltsC.sondeAlt./1000,30,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(compAltsC.meltAltInt./1000,compAltsC.sondeAlt./1000,30,'c','filled','MarkerEdgeColor','c');
scatter(compAltsC.meltAltMeas./1000,compAltsC.sondeAlt./1000,30,'b','filled','MarkerEdgeColor','b');
xlabel('HCR altitude (km)');
ylabel('Dropsonde altitude (km)');

text(minLim,maxLim+0.6*offset,'(a) CSET','fontweight','bold','fontsize',12);
grid on
%axis equal

xlim([minLim maxLim]);
ylim([minLim maxLim]);

measDiff=mean((compAltsC.meltAltMeas-compAltsC.sondeAlt)./1000,'omitnan');
measStd=std((compAltsC.meltAltMeas-compAltsC.sondeAlt)./1000,'omitnan');
intDiff=mean((compAltsC.meltAltInt-compAltsC.sondeAlt)./1000,'omitnan');
intStd=std((compAltsC.meltAltInt-compAltsC.sondeAlt)./1000,'omitnan');
estDiff=mean((compAltsC.meltAltEst-compAltsC.sondeAlt)./1000,'omitnan');
estStd=std((compAltsC.meltAltEst-compAltsC.sondeAlt)./1000,'omitnan');
zeroDiff=mean((compAltsC.zeroDegAlt-compAltsC.sondeAlt)./1000,'omitnan');
zeroStd=std((compAltsC.zeroDegAlt-compAltsC.sondeAlt)./1000,'omitnan');

text(minLim+offset/4,maxLim-0.5*offset,['Est.-drops.: ',num2str(estDiff,formatSpec), ...
    '\pm',num2str(estStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-1.5*offset,['0 deg.-drops.: ',num2str(zeroDiff,formatSpec), ...
    '\pm',num2str(zeroStd,formatSpec),' km'],'fontsize',12);

leg1=legend([l4,l3,l2,l1],{'Detections','Interpolations','Estimates','Zero deg'},'Location','southeast');
ax1.Position=[0.11 0.71 0.85 0.26];

ax2=subplot(3,1,2);
hold on;
maxLim=ceil(max(max(table2array(compAltsS(:,2:end))))/1000);
minLim=floor(min(min(table2array(compAltsS(:,2:end))))/1000);
offset=(maxLim-minLim)*0.08;
l0=plot([minLim maxLim],[minLim maxLim],'-r');
xticks([minLim:0.5:maxLim]);
yticks([minLim:0.5:maxLim]);

% Regression lines
fitOrth1=gmregress(compAltsS.zeroDegAlt./1000,compAltsS.sondeAlt./1000,1);
fitAll1=[fitOrth1(2) fitOrth1(1)];
xFitD =minLim:0.1:maxLim;
yFitD1 = polyval(fitAll1, xFitD);
l1=plot(xFitD, yFitD1,'-k','linewidth',1.5);
fitOrth2=gmregress(compAltsS.meltAltEst./1000,compAltsS.sondeAlt./1000,1);
fitAll2=[fitOrth2(2) fitOrth2(1)];
yFitD2 = polyval(fitAll2, xFitD);
l2=plot(xFitD, yFitD2,'color',[0.2 0.6 0.04],'linewidth',1.5);
fitOrth3=gmregress(compAltsS.meltAltInt./1000,compAltsS.sondeAlt./1000,1);
fitAll3=[fitOrth3(2) fitOrth3(1)];
yFitD3 = polyval(fitAll3, xFitD);
l3=plot(xFitD, yFitD3,'-c','linewidth',1.5);
fitOrth4=gmregress(compAltsS.meltAltMeas./1000,compAltsS.sondeAlt./1000,1);
fitAll4=[fitOrth4(2) fitOrth4(1)];
yFitD4 = polyval(fitAll4, xFitD);
l4=plot(xFitD, yFitD4,'-b','linewidth',1.5);

scatter(compAltsS.zeroDegAlt./1000,compAltsS.sondeAlt./1000,30,'k','filled','MarkerEdgeColor','k');
scatter(compAltsS.meltAltEst./1000,compAltsS.sondeAlt./1000,30,'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(compAltsS.meltAltInt./1000,compAltsS.sondeAlt./1000,30,'c','filled','MarkerEdgeColor','c');
scatter(compAltsS.meltAltMeas./1000,compAltsS.sondeAlt./1000,30,'b','filled','MarkerEdgeColor','b');
xlabel('HCR altitude (km)');
ylabel('Dropsonde altitude (km)');

text(minLim,maxLim+0.6*offset,'(b) SOCRATES','fontweight','bold','fontsize',12);
grid on
%axis equal

xlim([minLim maxLim]);
ylim([minLim maxLim]);

measDiff=mean((compAltsS.meltAltMeas-compAltsS.sondeAlt)./1000,'omitnan');
measStd=std((compAltsS.meltAltMeas-compAltsS.sondeAlt)./1000,'omitnan');
intDiff=mean((compAltsS.meltAltInt-compAltsS.sondeAlt)./1000,'omitnan');
intStd=std((compAltsS.meltAltInt-compAltsS.sondeAlt)./1000,'omitnan');
estDiff=mean((compAltsS.meltAltEst-compAltsS.sondeAlt)./1000,'omitnan');
estStd=std((compAltsS.meltAltEst-compAltsS.sondeAlt)./1000,'omitnan');
zeroDiff=mean((compAltsS.zeroDegAlt-compAltsS.sondeAlt)./1000,'omitnan');
zeroStd=std((compAltsS.zeroDegAlt-compAltsS.sondeAlt)./1000,'omitnan');

text(minLim+offset/4,maxLim-0.5*offset,['Det.-drops.: ',num2str(measDiff,formatSpec), ...
    '\pm',num2str(measStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-1.5*offset,['Int.-drops.: ',num2str(intDiff,formatSpec), ...
    '\pm',num2str(intStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-2.5*offset,['Est.-drops.: ',num2str(estDiff,formatSpec), ...
    '\pm',num2str(estStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-3.5*offset,['0 deg.-drops.: ',num2str(zeroDiff,formatSpec), ...
    '\pm',num2str(zeroStd,formatSpec),' km'],'fontsize',12);

leg1=legend([l4,l3,l2,l1],{'Detections','Interpolations','Estimates','Zero deg'},'Location','southeast');
ax2.Position=[0.11 0.385 0.85 0.26];

ax3=subplot(3,1,3);
hold on;
maxLim=ceil(max(max(table2array(compAltsO(:,2:end))))/1000);
minLim=floor(min(min(table2array(compAltsO(:,2:end))))/1000);
offset=(maxLim-minLim)*0.08;
l0=plot([minLim maxLim],[minLim maxLim],'-r');
xticks([minLim:0.5:maxLim]);
yticks([minLim:0.5:maxLim]);

% Regression lines
fitOrth1=gmregress(compAltsO.zeroDegAlt./1000,compAltsO.sondeAlt./1000,1);
fitAll1=[fitOrth1(2) fitOrth1(1)];
xFitD =minLim:0.1:maxLim;
yFitD1 = polyval(fitAll1, xFitD);
l1=plot(xFitD, yFitD1,'-k','linewidth',1.5);
fitOrth2=gmregress(compAltsO.meltAltEst./1000,compAltsO.sondeAlt./1000,1);
fitAll2=[fitOrth2(2) fitOrth2(1)];
yFitD2 = polyval(fitAll2, xFitD);
l2=plot(xFitD, yFitD2,'color',[0.2 0.6 0.04],'linewidth',1.5);
fitOrth3=gmregress(compAltsO.meltAltInt./1000,compAltsO.sondeAlt./1000,1);
fitAll3=[fitOrth3(2) fitOrth3(1)];
yFitD3 = polyval(fitAll3, xFitD);
l3=plot(xFitD, yFitD3,'-c','linewidth',1.5);
fitOrth4=gmregress(compAltsO.meltAltMeas./1000,compAltsO.sondeAlt./1000,1);
fitAll4=[fitOrth4(2) fitOrth4(1)];
yFitD4 = polyval(fitAll4, xFitD);
l4=plot(xFitD, yFitD4,'-b','linewidth',1.5);

scatter(compAltsO.zeroDegAlt./1000,compAltsO.sondeAlt./1000,30,'k','filled','MarkerEdgeColor','k');
scatter(compAltsO.meltAltEst./1000,compAltsO.sondeAlt./1000,30,[0.2 0.6 0.04],'MarkerEdgeColor',[0.2 0.6 0.04],'MarkerFaceColor',[0.2 0.6 0.04]);
scatter(compAltsO.meltAltInt./1000,compAltsO.sondeAlt./1000,30,'c','filled','MarkerEdgeColor','c');
scatter(compAltsO.meltAltMeas./1000,compAltsO.sondeAlt./1000,30,'b','filled','MarkerEdgeColor','b');
xlabel('HCR altitude (km)');
ylabel('Dropsonde altitude (km)');

text(minLim,maxLim+0.6*offset,'(c) OTREC','fontweight','bold','fontsize',12);
grid on
%axis equal

xlim([minLim maxLim]);
ylim([minLim maxLim]);

measDiff=mean((compAltsO.meltAltMeas-compAltsO.sondeAlt)./1000,'omitnan');
measStd=std((compAltsO.meltAltMeas-compAltsO.sondeAlt)./1000,'omitnan');
intDiff=mean((compAltsO.meltAltInt-compAltsO.sondeAlt)./1000,'omitnan');
intStd=std((compAltsO.meltAltInt-compAltsO.sondeAlt)./1000,'omitnan');
estDiff=mean((compAltsO.meltAltEst-compAltsO.sondeAlt)./1000,'omitnan');
estStd=std((compAltsO.meltAltEst-compAltsO.sondeAlt)./1000,'omitnan');
zeroDiff=mean((compAltsO.zeroDegAlt-compAltsO.sondeAlt)./1000,'omitnan');
zeroStd=std((compAltsO.zeroDegAlt-compAltsO.sondeAlt)./1000,'omitnan');

text(minLim+offset/4,maxLim-0.5*offset,['Det.-drops.: ',num2str(measDiff,formatSpec), ...
    '\pm',num2str(measStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-1.5*offset,['Int.-drops.: ',num2str(intDiff,formatSpec), ...
    '\pm',num2str(intStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-2.5*offset,['Est.-drops.: ',num2str(estDiff,formatSpec), ...
    '\pm',num2str(estStd,formatSpec),' km'],'fontsize',12);
text(minLim+offset/4,maxLim-3.5*offset,['0 deg.-drops.: ',num2str(zeroDiff,formatSpec), ...
    '\pm',num2str(zeroStd,formatSpec),' km'],'fontsize',12);

leg3=legend([l4,l3,l2,l1],{'Detections','Interpolations','Estimates','Zero deg'},'Location','southeast');
ax3.Position=[0.11 0.06 0.85 0.26];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'dropscatter'],'-dpng','-r0');
