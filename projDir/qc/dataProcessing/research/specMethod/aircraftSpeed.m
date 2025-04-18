% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

%% Speed bins
socrates=load('/scr/virga2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.2_full/spectralMoments/smoothingTest/socrates_smoothingAnalysis_everyOther_aircraftSpeed.mat');
otrec=load('/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/v3.2_full/spectralMoments/smoothingTest/otrec_smoothingAnalysis_everyOther_aircraftSpeed.mat');

%% Combine
socBoth=cat(2,socrates.edgesHalf',socrates.peakNum);
socBoth(any(isnan(socBoth),2),:)=[];
otBoth=cat(2,otrec.edgesHalf',otrec.peakNum);
otBoth(any(isnan(otBoth),2),:)=[];

comb=cat(1,socBoth,otBoth);
comb=sortrows(comb);


%% Fits

[pCs,sCs]=polyfit(socBoth(:,1),socBoth(:,2),2);
[yCs,dCs]=polyval(pCs,socBoth(:,1),sCs);

[pCo,sCo]=polyfit(otBoth(:,1),otBoth(:,2),2);
[yCo,dCo]=polyval(pCo,otBoth(:,1),sCo);

[pCb,sCb]=polyfit(comb(:,1),comb(:,2),2);
[yCb,dCb]=polyval(pCb,comb(:,1),sCb);

%% Scatter plot
close all

f1 = figure('Position',[200 500 800 400],'DefaultAxesFontSize',12,'renderer','painters');
t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

s1=nexttile(1);
hold on
l1=plot(socBoth(:,1),yCs,'-b','LineWidth',2);
l2=plot(socBoth(:,1),yCs+dCs,'--b','LineWidth',1);
plot(socBoth(:,1),yCs-dCs,'--b','LineWidth',1);
scatter(socBoth(:,1),socBoth(:,2),60,'^','filled','MarkerFaceColor','b','MarkerEdgeColor','w');
l3=plot(otBoth(:,1),yCo,'-g','LineWidth',2);
l4=plot(otBoth(:,1),yCo+dCo,'--g','LineWidth',1);
plot(otBoth(:,1),yCo-dCo,'--g','LineWidth',1);
scatter(otBoth(:,1),otBoth(:,2),60,'square','filled','MarkerFaceColor','g','MarkerEdgeColor','w');
xlim([5,275]);
ylim([4,21]);

legend([l1,l2,l3,l4],{, ...
    ['y=',num2str(pCs(1),3),'x^2-',num2str(abs(pCs(2)),3),'x+',num2str(pCs(3),3)], ...
    ['Mean standard error ',num2str(mean(dCs),2)], ...
    ['y=',num2str(pCb(1),3),'x^2-',num2str(abs(pCb(2)),3),'x+',num2str(pCb(3),3)], ...
    ['Mean standard error ',num2str(mean(dCb),2)]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Truncation point')

title(['(a) SOCRATES versus OTREC']);

grid on
box on

s2=nexttile(2);
hold on
l1=plot(comb(:,1),yCb,'-k','LineWidth',2);
l2=plot(comb(:,1),yCb+dCb,'--k','LineWidth',1);
plot(comb(:,1),yCb-dCb,'--k','LineWidth',1);
scatter(socBoth(:,1),socBoth(:,2),60,'^','filled','MarkerFaceColor','k','MarkerEdgeColor','w');
scatter(otBoth(:,1),otBoth(:,2),60,'square','filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([5,275]);
ylim([4,21]);

legend([l1,l2],{['y=',num2str(pCb(1),3),'x^2-',num2str(abs(pCb(2)),3),'x+',num2str(pCb(3),3)], ...
    ['Mean standard error ',num2str(mean(dCb),2)]});

xlabel('Aircraft speed (m s^{-1})')

title(['(b) SOCRATES and OTREC']);

grid on
box on

set(gcf,'PaperPositionMode','auto')
%print(f1,[figdir,'aircraftSpeed.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'aircraftSpeed.png'],'Resolution','300');