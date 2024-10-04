clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'oneError.mat']);
load([figdir,'aircSpeedHistSocrates.mat']);
load([figdir,'splitSig.mat']);

%% Figure

close all
numZero=2:250;
xlim1=[1,249];
xlim2=[400,800];

f1 = figure('Position',[200 500 600 1000],'DefaultAxesFontSize',12,'renderer','painters');

t = tiledlayout(3,1,'TileSpacing','tight','Padding','compact','TileIndexing', 'columnmajor');

s1=nexttile(1);

hold on

l0=plot(numZero,err21,'-b','LineWidth',2);
scatter(6,err21(5),'filled','MarkerFaceColor','red','MarkerEdgeColor','k');
xlabel('Truncation piont');
ylabel('RMSE (dB)');
ylim([5,11]);
xlim(xlim1);

grid on
box on
title('(a) Root mean square error of example spectrum')
annotation('textarrow',[0.27,0.13],[0.84,0.75],'String','Optimum truncation point','FontSize',12)

s2=nexttile(2);

bar(numZero(1:end-1),H./sum(H).*100,1,'red')

xlabel('Truncation piont');
ylabel('Frequency (%)')
title('(b) Distribution of optimum truncation points')
annotation('textarrow',[0.27,0.135],[0.615,0.6],'String',' Most frequent optimum truncation point','FontSize',12)

grid on
box on

s3=nexttile(3);

hold on

l1=plot(1:2:length(split.filt7)-1,split.sig1,'-k','LineWidth',1);
l2=plot(1:2:length(split.filt7)-1,split.sig2,'-','LineWidth',1,'Color',[0.6,0.6,0.6]);
l3=plot(1:2:length(split.filt7)-1,split.filt40,'-','LineWidth',2,'Color',[0.8,0,0.8]);
l4=plot(1:length(split.filt7),split.filt7,'-g','LineWidth',2);
ylabel('Power (dB)');
xlabel('Spectral bin number')
ylim([-75,-25]);
xlim(xlim2);

grid on
box on
title('(c) Example with optimum truncation point of 40')
legend([l1,l2,l3,l4],{'Signal A';'Signal B';'Filtered TP=40';'Filtered TP=7'},'Location','northwest');

h1=drawellipse('Center',[530,-53],'SemiAxes',[15,9.5],'Color','b', ...
    'FaceAlpha',0,'HandleVisibility','off','MarkerSize',0.01);
h2=drawellipse('Center',[740,-50],'SemiAxes',[12,8],'Color','b', ...
    'FaceAlpha',0,'HandleVisibility','off','MarkerSize',0.01);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'optTP.png'],'-dpng','-r0');
