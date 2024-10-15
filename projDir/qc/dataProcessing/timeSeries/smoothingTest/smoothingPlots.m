% Analyze HCR time series
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/spectralMoments/smoothingTest/projectComparison/';

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

[pLs,sLs]=polyfit(socBoth(:,1),socBoth(:,2),1);
[yLs,dLs]=polyval(pLs,socBoth(:,1),sLs);

[pCs,sCs]=polyfit(socBoth(:,1),socBoth(:,2),2);
[yCs,dCs]=polyval(pCs,socBoth(:,1),sCs);

[pLo,sLo]=polyfit(otBoth(:,1),otBoth(:,2),1);
[yLo,dLo]=polyval(pLo,otBoth(:,1),sLo);

[pCo,sCo]=polyfit(otBoth(:,1),otBoth(:,2),2);
[yCo,dCo]=polyval(pCo,otBoth(:,1),sCo);

[pLb,sLb]=polyfit(comb(:,1),comb(:,2),1);
[yLb,dLb]=polyval(pLb,comb(:,1),sLb);

[pCb,sCb]=polyfit(comb(:,1),comb(:,2),2);
[yCb,dCb]=polyval(pCb,comb(:,1),sCb);

%% Scatter plot
close all

f1 = figure('Position',[200 500 1200 800],'DefaultAxesFontSize',12,'renderer','painters');
t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
l1=plot(socBoth(:,1),yLs,'-r','LineWidth',2);
l2=plot(socBoth(:,1),yLs+dLs,'--r','LineWidth',1);
plot(socBoth(:,1),yLs-dLs,'--r','LineWidth',1);
l3=plot(socBoth(:,1),yCs,'-g','LineWidth',2);
l4=plot(socBoth(:,1),yCs+dCs,'--g','LineWidth',1);
plot(socBoth(:,1),yCs-dCs,'--g','LineWidth',1);
scatter(socBoth(:,1),socBoth(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(socBoth(:,1))-1,max(socBoth(:,1))+1]);
ylim([min(socBoth(:,2))-1,max(socBoth(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLs(1)),'x+',num2str(pLs(2))], ...
    ['Mean standard error ',num2str(mean(dLs))], ...
    ['y=',num2str(pCs(1)),'x^2+',num2str(pCs(2)),'x+',num2str(pCs(3))], ...
    ['Mean standard error ',num2str(mean(dCs))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['SOCRATES binned']);

grid on
box on

s2=nexttile(2);
hold on
l1=plot(otBoth(:,1),yLo,'-r','LineWidth',2);
l2=plot(otBoth(:,1),yLo+dLo,'--r','LineWidth',1);
plot(otBoth(:,1),yLo-dLo,'--r','LineWidth',1);
l3=plot(otBoth(:,1),yCo,'-g','LineWidth',2);
l4=plot(otBoth(:,1),yCo+dCo,'--g','LineWidth',1);
plot(otBoth(:,1),yCo-dCo,'--g','LineWidth',1);
scatter(otBoth(:,1),otBoth(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(otBoth(:,1))-1,max(otBoth(:,1))+1]);
ylim([min(otBoth(:,2))-1,max(otBoth(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLo(1)),'x+',num2str(pLo(2))], ...
    ['Mean standard error ',num2str(mean(dLo))], ...
    ['y=',num2str(pCo(1)),'x^2+',num2str(pCo(2)),'x+',num2str(pCo(3))], ...
    ['Mean standard error ',num2str(mean(dCo))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['OTREC binned']);

grid on
box on

s3=nexttile(3);
hold on
l1=plot(comb(:,1),yLb,'-r','LineWidth',2);
l2=plot(comb(:,1),yLb+dLb,'--r','LineWidth',1);
plot(comb(:,1),yLb-dLb,'--r','LineWidth',1);
l3=plot(comb(:,1),yCb,'-g','LineWidth',2);
l4=plot(comb(:,1),yCb+dCb,'--g','LineWidth',1);
plot(comb(:,1),yCb-dCb,'--g','LineWidth',1);
scatter(comb(:,1),comb(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(comb(:,1))-1,max(comb(:,1))+1]);
ylim([min(comb(:,2))-1,max(comb(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLb(1)),'x+',num2str(pLb(2))], ...
    ['Mean standard error ',num2str(mean(dLb))], ...
    ['y=',num2str(pCb(1)),'x^2+',num2str(pCb(2)),'x+',num2str(pCb(3))], ...
    ['Mean standard error ',num2str(mean(dCb))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['SOCRATES and OTREC binned']);

grid on
box on

%% Cases
socrates=load('/scr/virga2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.2_full/spectralMoments/smoothingTest/socrates_smoothingAnalysis_everyOther_aircraftSpeed_cases.mat');
otrec=load('/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/v3.2_full/spectralMoments/smoothingTest/otrec_smoothingAnalysis_everyOther_aircraftSpeed_cases.mat');

%% Combine
socBoth=cat(2,socrates.velAircMeanC,socrates.peakNumCases);
socBoth(any(isnan(socBoth),2),:)=[];
otBoth=cat(2,otrec.velAircMeanC,otrec.peakNumCases);
otBoth(any(isnan(otBoth),2),:)=[];

comb=cat(1,socBoth,otBoth);
comb=sortrows(comb);


%% Fits

[pLs,sLs]=polyfit(socBoth(:,1),socBoth(:,2),1);
[yLs,dLs]=polyval(pLs,socBoth(:,1),sLs);

[pCs,sCs]=polyfit(socBoth(:,1),socBoth(:,2),2);
[yCs,dCs]=polyval(pCs,socBoth(:,1),sCs);

[pLo,sLo]=polyfit(otBoth(:,1),otBoth(:,2),1);
[yLo,dLo]=polyval(pLo,otBoth(:,1),sLo);

[pCo,sCo]=polyfit(otBoth(:,1),otBoth(:,2),2);
[yCo,dCo]=polyval(pCo,otBoth(:,1),sCo);

[pLb,sLb]=polyfit(comb(:,1),comb(:,2),1);
[yLb,dLb]=polyval(pLb,comb(:,1),sLb);

[pCb,sCb]=polyfit(comb(:,1),comb(:,2),2);
[yCb,dCb]=polyval(pCb,comb(:,1),sCb);

%% Scatter plot

s4=nexttile(4);
hold on
l1=plot(socBoth(:,1),yLs,'-r','LineWidth',2);
l2=plot(socBoth(:,1),yLs+dLs,'--r','LineWidth',1);
plot(socBoth(:,1),yLs-dLs,'--r','LineWidth',1);
l3=plot(socBoth(:,1),yCs,'-g','LineWidth',2);
l4=plot(socBoth(:,1),yCs+dCs,'--g','LineWidth',1);
plot(socBoth(:,1),yCs-dCs,'--g','LineWidth',1);
scatter(socBoth(:,1),socBoth(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(socBoth(:,1))-1,max(socBoth(:,1))+1]);
ylim([min(socBoth(:,2))-1,max(socBoth(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLs(1)),'x+',num2str(pLs(2))], ...
    ['Mean standard error ',num2str(mean(dLs))], ...
    ['y=',num2str(pCs(1)),'x^2+',num2str(pCs(2)),'x+',num2str(pCs(3))], ...
    ['Mean standard error ',num2str(mean(dCs))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['SOCRATES cases']);

grid on
box on

s5=nexttile(5);
hold on
l1=plot(otBoth(:,1),yLo,'-r','LineWidth',2);
l2=plot(otBoth(:,1),yLo+dLo,'--r','LineWidth',1);
plot(otBoth(:,1),yLo-dLo,'--r','LineWidth',1);
l3=plot(otBoth(:,1),yCo,'-g','LineWidth',2);
l4=plot(otBoth(:,1),yCo+dCo,'--g','LineWidth',1);
plot(otBoth(:,1),yCo-dCo,'--g','LineWidth',1);
scatter(otBoth(:,1),otBoth(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(otBoth(:,1))-1,max(otBoth(:,1))+1]);
ylim([min(otBoth(:,2))-1,max(otBoth(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLo(1)),'x+',num2str(pLo(2))], ...
    ['Mean standard error ',num2str(mean(dLo))], ...
    ['y=',num2str(pCo(1)),'x^2+',num2str(pCo(2)),'x+',num2str(pCo(3))], ...
    ['Mean standard error ',num2str(mean(dCo))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['OTREC cases']);

grid on
box on

s6=nexttile(6);
hold on
l1=plot(comb(:,1),yLb,'-r','LineWidth',2);
l2=plot(comb(:,1),yLb+dLb,'--r','LineWidth',1);
plot(comb(:,1),yLb-dLb,'--r','LineWidth',1);
l3=plot(comb(:,1),yCb,'-g','LineWidth',2);
l4=plot(comb(:,1),yCb+dCb,'--g','LineWidth',1);
plot(comb(:,1),yCb-dCb,'--g','LineWidth',1);
scatter(comb(:,1),comb(:,2),60,'filled','MarkerFaceColor','k','MarkerEdgeColor','w');
xlim([min(comb(:,1))-1,max(comb(:,1))+1]);
ylim([min(comb(:,2))-1,max(comb(:,2))+6]);

legend([l1,l2,l3,l4],{['y=',num2str(pLb(1)),'x+',num2str(pLb(2))], ...
    ['Mean standard error ',num2str(mean(dLb))], ...
    ['y=',num2str(pCb(1)),'x^2+',num2str(pCb(2)),'x+',num2str(pCb(3))], ...
    ['Mean standard error ',num2str(mean(dCb))]});

xlabel('Aircraft speed (m s^{-1})')
ylabel('Peak at number of non-zeros')

title(['SOCRATES and OTREC cases']);

grid on
box on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'smoothingAnalysis_everyOther_aircraftSpeed'],'-dpng','-r0');
