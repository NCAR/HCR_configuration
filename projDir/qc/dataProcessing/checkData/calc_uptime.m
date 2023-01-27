% Compare different data sets to make sure the data is complete

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
qcGood='qc1'; % field, qc0, qc1, qc2
qcVersionGood='v1.1';
freqGood='10hz'; % 10hz, 2hz, 2hzMerged

thresh12=[10 50];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_uptime.txt'];

caseList = table2array(readtable(infile));

indirGood=HCRdir(project,qcGood,qcVersionGood,freqGood);

figdir=[indirGood(1:end-5),'uptime/'];

%% Run processing

shouldTime=[];
isTime=[];

% Go through flights
for ii=1:size(caseList,1)
    disp(['Flight ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    %endTime=startTime+hours(1);

    fileListGood=makeFileList(indirGood,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    %% Load data

    data=[];
    data.FLAG=[];

    % Load data
    data=read_HCR(fileListGood,data,startTime,endTime);


    % Resample

    startTimeReal=dateshift(data.time(1),'start','second');
    endTimeReal=dateshift(data.time(end),'end','second');

    newTime=startTimeReal:seconds(0.1):endTimeReal;

    missingInds=any(data.FLAG==11,1);

    shouldTime=cat(1,shouldTime,length(newTime));

    isTime=cat(1,isTime,length(data.time)-sum(missingInds));
end

%% Calc vars
upPerc=isTime./shouldTime.*100;
upPercTot=sum(isTime)/sum(shouldTime)*100;

shouldHours=shouldTime./10./60./60;
isHours=isTime./10./60./60;

shouldIs=cat(2,shouldHours,isHours);
if size(shouldIs,1)==1
    shouldIs=cat(1,shouldIs,[nan,nan]);
end
%% Plot

close all

f1=figure('DefaultAxesFontSize',12,'renderer','painters');
set(f1,'Position',[0.5 0.1 800 600]);

s1=subplot(2,1,1);
bar(shouldIs);
xlim([0.5,length(shouldHours)+0.5]);
xticks(1:length(shouldHours));
ylims=s1.YLim;
ylim([0,ylims(2)+2]);
legend('Requested','Operated','Location','Northeast','Orientation','horizontal');
set(gca, 'YGrid', 'on', 'XGrid', 'off');
set(gca,'TickLength',[0,0]);

xlabel('Flight number')
ylabel('Hours')

title([project,' HCR operation hours']);

s2=subplot(2,1,2);
bar(upPerc);
xlim([0.5,length(shouldHours)+0.5]);
xticks(1:length(shouldHours));
ylim([0,100]);
set(gca, 'YGrid', 'on', 'XGrid', 'off');
set(gca,'TickLength',[0,0]);

xlabel('Flight number')
ylabel('Uptime (%)')

title(['HCR uptime. Total uptime: ',num2str(upPercTot,3)]);
set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_operationHours.png'],'-dpng','-r0')

save([figdir,project,'_operationHours.mat'],'isHours','shouldHours');
