clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'specParamsDemo.mat']);

%% Prepare

powerAbove=demSpec.powerSmooth;
powerAbove(demSpec.powerSmooth<demSpec.noiseF)=nan;

leftEdge=find(~isnan(powerAbove),1,'first');
rightEdge=find(~isnan(powerAbove),1,'last');

peaks=islocalmax(powerAbove);
peaks=find(peaks==1);
leftPeak=peaks(1);
rightPeak=peaks(2);

leftEdgeI=[demSpec.vel(leftEdge),demSpec.noiseF(leftEdge)];
rightEdgeI=[demSpec.vel(rightEdge),demSpec.noiseF(rightEdge)];
leftPeakI=[demSpec.vel(leftPeak),demSpec.powerSmooth(leftPeak)];
rightPeakI=[demSpec.vel(rightPeak),demSpec.powerSmooth(rightPeak)];

e2eHalf=(rightEdge-leftEdge)/2;
e2eMid=leftEdge+e2eHalf;
halfDist=(rightEdgeI(1)-leftEdgeI(1))/2;

%% Figure

close all

f1 = figure('Position',[200 500 600 350],'DefaultAxesFontSize',12,'renderer','painters');

t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact','TileIndexing', 'columnmajor');

s1=nexttile(1);

hold on

plot(demSpec.vel,demSpec.powerSmooth,'-','LineWidth',1.5,'Color',[0.7,0.7,0.7]);
l1=plot(demSpec.vel,powerAbove,'-k','LineWidth',2);
l2=plot(demSpec.vel,demSpec.noiseF,'-c','LineWidth',1.5);

plot([leftEdgeI(1),leftPeakI(1)],[leftEdgeI(2),leftPeakI(2)],'-','Color',[1 0 1],'LineWidth',1.5);
plot([leftPeakI(1),rightEdgeI(1)],[leftPeakI(2),rightEdgeI(2)],'-','Color',[0.6940 0.3840 0.7560],'LineWidth',1.5);

e=errorbar(demSpec.vel(e2eMid),demSpec.noiseF(1)-4,halfDist,'horizontal');
e.Color=[0.4660 0.6740 0.1880];
e.LineWidth=1.5;

scatter(leftEdgeI(1),leftEdgeI(2),60,'filled','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerEdgeColor','k');
scatter(rightEdgeI(1),rightEdgeI(2),60,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor','k');
scatter(leftPeakI(1),leftPeakI(2),60,'filled','MarkerFaceColor',[0 0.6470 0.9410],'MarkerEdgeColor','k');
scatter(rightPeakI(1),rightPeakI(2),60,'filled','MarkerFaceColor',[0.4010 0.8450 0.7],'MarkerEdgeColor','k');

xlabel('Velocity (m s^{-1})');
ylabel('Power (dB)');
ylim([-60,10]);
xlim([demSpec.vel(1),demSpec.vel(end)]);

grid on
box on
legend([l1,l2],{'Power spectrum';'Spec. noise floor'},'Location','northeast');

t = annotation('textbox',[0.14, 0.3, 0.1, 0.1],'String',"Left edge",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0.9290 0.6940 0.1250]);
t = annotation('textbox',[0.80, 0.3, 0.1, 0.1],'String',"Right edge",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0.8500 0.3250 0.0980]);
t = annotation('textbox',[0.31, 0.8, 0.1, 0.1],'String',"Left peak",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0 0.6470 0.9410]);
t = annotation('textbox',[0.67, 0.7, 0.1, 0.1],'String',"Right peak",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0.4010 0.8450 0.7]);
t = annotation('textbox',[0.41, 0.4, 0.1, 0.1],'String',"Left slope",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[1 0 1],'Rotation',60);
t = annotation('textbox',[0.54, 0.58, 0.1, 0.1],'String',"Right slope",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0.6940 0.3840 0.7560],'Rotation',-44);
t = annotation('textbox',[0.46, 0.16, 0.1, 0.1],'String',"Edge-to-edge width",'FontSize',12,'Margin',1, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','BackgroundColor',[0.4660 0.6740 0.1880]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'specParamsDemo.png'],'-dpng','-r0');
