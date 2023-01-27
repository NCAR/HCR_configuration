% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figDirs={'/scr/snow2/rsfdata/projects/noreaster/hcr/qc2/cfradial/v2.0_full/uptime/';
    '/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/uptime/';
    '/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.1_full/uptime/';
    '/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/v3.1_full/uptime/';
    '/scr/sleet3/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.1_full/uptime/'};

projects={'noreaster';'cset';'socrates';'otrec';'spicule'};

figdir=figDirs{1};

%% Run processing

shouldHoursAll=[];
isHoursAll=[];
upPercProj=[];
totHoursAll=[];

% Go through flights
for ii=1:length(projects)
    loadIn=load([figDirs{ii},projects{ii},'_operationHours.mat']);
    
    upPercProj=cat(1,upPercProj,sum(loadIn.isHours)/sum(loadIn.shouldHours)*100);

    shouldHoursAll=cat(1,shouldHoursAll,loadIn.shouldHours);
    isHoursAll=cat(1,isHoursAll,loadIn.isHours);

    totHoursAll=cat(1,totHoursAll,sum(loadIn.isHours));

end

%% Calc vars
upPerc=isHoursAll./shouldHoursAll.*100;
upPercTot=sum(isHoursAll)/sum(shouldHoursAll)*100;

%% Plot

close all

f1=figure('DefaultAxesFontSize',12,'renderer','painters');
set(f1,'Position',[0.5 0.1 600 400]);

bar(upPercProj);
xlim([0.5,length(upPercProj)+0.5]);
xticks(1:length(upPercProj));
xticklabels(projects);
ylim([0,100]);
set(gca, 'YGrid', 'on', 'XGrid', 'off');
set(gca,'TickLength',[0,0]);

ylabel('Uptime (%)')

text(1:length(upPercProj),repmat(10,1,length(upPercProj)),num2str(totHoursAll,4),...
    'FontSize',12,'HorizontalAlignment','center','BackgroundColor','w');

text(1:length(upPercProj),repmat(90,1,length(upPercProj)),num2str(upPercProj,3),...
    'FontSize',12,'HorizontalAlignment','center','BackgroundColor','w');

title(['Total uptime: ',num2str(upPercTot,3),'. Total hours: ',num2str(sum(totHoursAll),4),'.']);
set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'operationHours.png'],'-dpng','-r0')
