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

figdir=[dataDir(1:end-5),'cloudClass/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/cloudClass_',project,'.txt'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

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

    %% Un-truncate

    disp('Adding gaps back in ...')
    data=unJoinOverMissing(data,nonMissingInds);
       
    %% Cloud classification

    disp('Cloud classification ...')
    cloudClass=findCloudClass(data.ECHO_TYPE_2D,data.cloudPuzzle_echoType,data.TEMP,data.elevation,data.TOPO,data.asl);

    % Fill in small with not classified
    cloudClass(~isnan(data.ECHO_TYPE_2D) & isnan(cloudClass))=0;
       
    %% Plot

    indWant=3000;
    indSpacing=round(length(data.time)/indWant);
    getInds=1:indSpacing:length(data.time);

    plotTime=data.time(getInds);
    plotAsl=data.asl(:,getInds);

    % Prepare cloud class
    % Not classified=0
    % stratLow=1
    % stratMid=2;
    % stratHigh=3
    % stratPrecipShallow=11
    % stratPrecipMid=12
    % stratPrecipDeep=13
    % convYoungShallow=21
    % convYoungMid=22
    % convYongDeep=23
    % convMatureShallow=31
    % convMatureMid=32
    % convMatureDeep=33

    classPlot=cloudClass(:,getInds);
    classPlot(classPlot==11)=4;
    classPlot(classPlot==12)=5;
    classPlot(classPlot==13)=6;
    classPlot(classPlot==21)=7;
    classPlot(classPlot==22)=8;
    classPlot(classPlot==23)=9;
    classPlot(classPlot==31)=10;
    classPlot(classPlot==32)=11;
    classPlot(classPlot==33)=12;

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

    % Prepare strat conv
    
    disp('Plotting ...');
    
    close all
    
    csSubPlot=data.ECHO_TYPE_2D(:,getInds);
    csSubPlot(csSubPlot==14)=1;
    csSubPlot(csSubPlot==16)=2;
    csSubPlot(csSubPlot==18)=3;
    csSubPlot(csSubPlot==25)=4;
    csSubPlot(csSubPlot==30)=5;
    csSubPlot(csSubPlot==32)=6;
    csSubPlot(csSubPlot==34)=7;
    csSubPlot(csSubPlot==36)=8;
    csSubPlot(csSubPlot==38)=9;

    colmapSC=[0,0.1,0.6;
        0.38,0.42,0.96;
        0.65,0.74,0.86;
        0.32,0.78,0.59;
        0.7,0,0;
        1,0,1;
        1,1,0;
        0.99,0.77,0.22;
        1,0,0];

    close all
    
    f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible',showPlot);
       
    s1=subplot(2,1,1);
    
    hold on
    surf(plotTime,plotAsl./1000,csSubPlot,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    s1.Colormap=colmapSC;
    caxis([0.5 9.5]);
    cb=colorbar;
    cb.Ticks=1:9;
    cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
        'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
    grid on
    box on
    title('Echo type')

    s2=subplot(2,1,2);
    
    hold on
    surf(plotTime,plotAsl./1000,classPlot,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-1 13]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    s2.Colormap=colmapCC;
    caxis([-0.5 12.5]);
    cb=colorbar;
    cb.Ticks=0:12;
    cb.TickLabels={'Not classified','Strat Low','Strat Mid','Strat High',...
        'Strat Precip Shallow','Strat Precip Mid','Strat Precip Deep',...
        'Conv Young Shallow','Conv Young Mid','Conv Yong Deep',...
        'Conv Mature Shallow','Conv Mature Mid','Conv Mature Deep'};
    grid on
    box on
    title('Cloud classification')
        
    linkaxes([s1],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_cloudClass_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
   
end