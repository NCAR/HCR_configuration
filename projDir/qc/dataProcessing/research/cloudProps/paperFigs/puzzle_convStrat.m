% Call cloud classification script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

showPlot='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

% Loop through cases

caseList=[2019,9,4,12,36,2019,9,4,12,52; ...
    2019,9,22,14,28,2019,9,22,14,38; ...
    2019,8,18,17,18,2019,8,18,18,11];
caseStart=datetime([caseList(:,1:5),zeros(3,1)]);
caseEnd=datetime([caseList(:,6:10),zeros(3,1)]);

puzz.c1=[];
puzz.c2=[];
puzz.c3=[];

timeC.c1=[];
timeC.c2=[];
timeC.c3=[];

aslC.c1=[];
aslC.c2=[];
aslC.c3=[];

cloudCountAll=[];

fieldsP=fields(puzz);

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Get data

    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.ECHO_TYPE_2D=[];
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    ylimUpper=(max(data.asl(~isnan(data.ECHO_TYPE_2D)))./1000)+0.5;

    echoPlot=nan(size(data.ECHO_TYPE_2D));
    echoPlot(data.ECHO_TYPE_2D<20)=1;
    echoPlot(data.ECHO_TYPE_2D>20 & data.ECHO_TYPE_2D<30)=2;
    echoPlot(data.ECHO_TYPE_2D>=30)=3;

    %% Truncate to non missing

    puzz.(fieldsP{aa})=echoPlot;
    timeC.(fieldsP{aa})=data.time;
    aslC.(fieldsP{aa})=data.asl;
end

%% Plot

disp('Plotting ...');

close all

fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1200,600]);

ax1=subplot(2,2,1);

% Get indices
indWant=3000;
indSpacing=round(length(timeC.c1)/indWant);
getInds=1:indSpacing:length(timeC.c1);

plotTime=timeC.c1(getInds);
plotASL=aslC.c1(:,getInds);
plotPuzzle=puzz.c1(:,getInds);


colMap=[0,0,1;
    0.32,0.78,0.59;
    0.7,0,0];

hold on;
sub2=surf(plotTime,plotASL./1000,plotPuzzle,'edgecolor','none');
view(2);
ax1.Colormap=colMap;
caxis([0.5 3.5]);
ylim([0,ylimUpper]);
ylabel('Altitude (km)');
xlim([plotTime(1),plotTime(end)]);
%title('Cloud Puzzle')
grid on
box on

%text(datetime(2019,9,4,12,36,15),13.4,'(a)','FontSize',14,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

ax2=subplot(2,2,2);

% Get indices
indWant=3000;
indSpacing=round(length(timeC.c2)/indWant);
getInds=1:indSpacing:length(timeC.c2);

plotTime=timeC.c2(getInds);
plotASL=aslC.c2(:,getInds);
plotPuzzle=puzz.c2(:,getInds);

hold on;
sub2=surf(plotTime,plotASL./1000,plotPuzzle,'edgecolor','none');
view(2);
ax2.Colormap=colMap;
caxis([0.5 3.5]);
cb=colorbar;
cb.Ticks=1:3;
cb.TickLabels={'Stratiform','Mixed','Convective'};
ylim([0,ylimUpper]);
xlim([plotTime(1),plotTime(end)]);
%title('Cloud Puzzle')
grid on
box on

%text(datetime(2019,9,22,14,28,15),13.4,'(b)','FontSize',14,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

ax2.SortMethod='childorder';

ax3=subplot(2,1,2);

% Get indices
indWant=3000;
indSpacing=round(length(timeC.c3)/indWant);
getInds=1:indSpacing:length(timeC.c3);

plotTime=timeC.c3(getInds);
plotASL=aslC.c3(:,getInds);
plotPuzzle=puzz.c3(:,getInds);

hold on;
sub3=surf(plotTime,plotASL./1000,plotPuzzle,'edgecolor','none');
view(2);
ax3.Colormap=colMap;
caxis([0.5 3.5]);
ylim([0,ylimUpper]);
ylabel('Altitude (km)');
xlim([plotTime(1),plotTime(end)]);
%title('Cloud Puzzle')
grid on
box on

%text(datetime(2019,8,18,17,18,30),13.4,'(c)','FontSize',14,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

ax3.SortMethod='childorder';

ax1.Position=[0.047,0.57,0.4,0.41];
ax2.Position=[0.486,0.57,0.42,0.41];
ax3.Position=[0.047,0.08,0.92,0.41];
cb.Position=[0.91,0.57,0.0178,0.41];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'cloudPuzzle_convStrat.png'],'-dpng','-r0');
