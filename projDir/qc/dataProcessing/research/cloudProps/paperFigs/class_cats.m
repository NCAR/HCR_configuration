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

caseList=[2019,8,12,12,30,2019,8,12,12,50; ...
    2019,8,16,15,15,2019,8,16,15,54];
caseStart=datetime([caseList(:,1:5),zeros(2,1)]);
caseEnd=datetime([caseList(:,6:10),zeros(2,1)]);

classC.c1=[];
classC.c2=[];

timeC.c1=[];
timeC.c2=[];

aslC.c1=[];
aslC.c2=[];

cloudCountAll=[];

fieldsP=fields(classC);

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Get data

    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ=[];
    data.ECHO_TYPE_2D=[];
    data.FLAG=[];
    data.ANTFLAG=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    data.TOPO=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    ylimUpper=(max(data.asl(~isnan(data.ECHO_TYPE_2D)))./1000)+0.5;

    %% Truncate to non missing

    gapSecs=10;
    [data,nonMissingInds]=joinOverMissing(data,gapSecs);

    %% Create cloudID

    minCloudSizePix=1000;

    cloudID=makeCloudID(data,minCloudSizePix);

    %% Cloud puzzle
    % Breaks up really big clouds

    disp('Breaking up large ...')
    data.cloudPuzzle=f_cloudPuzzle_breakLarge(cloudID,data);

    % Breaks out isolated convective that penetrate into stratiform
    disp('Breaking out isolated conv ...')
    data.cloudPuzzle_echoType=breakout_isolatedConv(data);

    uClouds=unique(data.cloudPuzzle(~isnan(data.cloudPuzzle)));
    cloudCount=length(uClouds);

    cloudCountAll=[cloudCountAll,cloudCount];

    %% Un-truncate

    disp('Adding gaps back in ...')
    data=unJoinOverMissing(data,nonMissingInds);

     %% Cloud classification

    disp('Cloud classification ...')
    cloudClass=findCloudClass(data.ECHO_TYPE_2D,data.cloudPuzzle_echoType,data.TEMP,data.MELTING_LAYER,data.elevation,data.TOPO,data.asl);

    % Fill in small with not classified
    cloudClass(~isnan(data.ECHO_TYPE_2D) & isnan(cloudClass))=0;

    classC.(fieldsP{aa})=cloudClass;
    timeC.(fieldsP{aa})=data.time;
    aslC.(fieldsP{aa})=data.asl;
end

%% Plot

colmapCC=[0,0,0;
    204,255,204;
    153,204,0;
    0,128,0;
    0,204,255;
    51,102,255;
    0,0,180;
    255,204,0;
    255,102,0;
    220,0,0;
    255,153,220;
    204,153,255;
    128,0,128];

colmapCC=colmapCC./255;

disp('Plotting ...');

close all

fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1200,600]);

ax1=subplot(2,1,1);

% Get indices
indWant=3000;
indSpacing=round(length(timeC.c1)/indWant);
getInds=1:indSpacing:length(timeC.c1);

plotTime=timeC.c1(getInds);
plotASL=aslC.c1(:,getInds);
classPlot=classC.c1(:,getInds);

classPlot(classPlot==11)=4;
classPlot(classPlot==12)=5;
classPlot(classPlot==13)=6;
classPlot(classPlot==21)=7;
classPlot(classPlot==22)=8;
classPlot(classPlot==23)=9;
classPlot(classPlot==31)=10;
classPlot(classPlot==32)=11;
classPlot(classPlot==33)=12;

hold on;
sub2=surf(plotTime,plotASL./1000,classPlot,'edgecolor','none');
view(2);
ax1.Colormap=colmapCC;
caxis([-0.5 12.5]);
ylim([0,11.5]);
ylabel('Altitude (km)');
xlim([plotTime(1),plotTime(end)]);
%title('Cloud Puzzle')
grid on
box on

text(datetime(2019,8,12,12,30,15),10.4,'(a)','FontSize',14,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

ax2=subplot(2,1,2);

% Get indices
indWant=3000;
indSpacing=round(length(timeC.c2)/indWant);
getInds=1:indSpacing:length(timeC.c2);

plotTime=timeC.c2(getInds);
plotASL=aslC.c2(:,getInds);
classPlot=classC.c2(:,getInds);

classPlot(classPlot==11)=4;
classPlot(classPlot==12)=5;
classPlot(classPlot==13)=6;
classPlot(classPlot==21)=7;
classPlot(classPlot==22)=8;
classPlot(classPlot==23)=9;
classPlot(classPlot==31)=10;
classPlot(classPlot==32)=11;
classPlot(classPlot==33)=12;

hold on;
sub2=surf(plotTime,plotASL./1000,classPlot,'edgecolor','none');
view(2);
ax2.Colormap=colmapCC;
caxis([-0.5 12.5]);
cb=colorbar;
cb.Ticks=0:12;
cb.TickLabels={'Not classified','CloudLow','CloudMid','CloudHigh',...
    'StratShallow','StratMid','StratDeep',...
    'ConvShallow','ConvMid','ConvDeep',...
    'ConvStratShallow','ConvStratMid','ConvStratDeep'};
ylim([0,14]);
xlim([plotTime(1),plotTime(end)]);
ylabel('Altitude (km)');
%title('Cloud Puzzle')
grid on
box on

text(datetime(2019,8,16,15,15,30),12.6,'(b)','FontSize',14,'FontWeight','bold','BackgroundColor','w','EdgeColor','k');

ax2.SortMethod='childorder';

ax1.Position=[0.047,0.57,0.81,0.41];
ax2.Position=[0.047,0.08,0.81,0.41];

cb.Position=[0.865,0.1,0.018,0.86];

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'cloudClassCats.png'],'-dpng','-r0');
