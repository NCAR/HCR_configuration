% Call cloud puzzle script

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

figdir=[dataDir(1:end-5),'cloudPuzzleEchoType/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/cloudPuzzle_',project,'.txt'];

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
        
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
 
    ylimUpper=(max(data.asl(data.FLAG==1))./1000)+0.5;
   
    %% Truncate to non missing

    gapSecs=10;
    [data,nonMissingInds]=joinOverMissing(data,gapSecs);

    %% Create cloudID

    minCloudSizePix=1000;

    cloudID=makeCloudID(data,minCloudSizePix);

    %% Cloud puzzle
    % Breaks up really big clouds
    % Breaks out convective regions that insert into stratiform regions

    data.cloudPuzzle=cloudID;
    data.cloudPuzzle(data.cloudPuzzle==0)=nan;

    uClouds=unique(data.cloudPuzzle(~isnan(data.cloudPuzzle)));
    cloudCount=length(uClouds);

    %% Un-truncate

    data=unJoinOverMissing(data,nonMissingInds);
       
    %% Plot
    
    disp('Plotting ...');

    % Get indices
    indWant=3000;
    indSpacing=round(length(data.time)/indWant);
    getInds=1:indSpacing:length(data.time);

    plotTime=data.time(getInds);
    plotASL=data.asl(:,getInds);
    plotDBZ=data.DBZ(:,getInds);
    plotPuzzle=data.cloudPuzzle(:,getInds);

    csPlot=data.ECHO_TYPE_2D(:,getInds);
    csPlot(csPlot==14)=1;
    csPlot(csPlot==16)=2;
    csPlot(csPlot==18)=3;
    csPlot(csPlot==25)=4;
    csPlot(csPlot==30)=5;
    csPlot(csPlot==32)=6;
    csPlot(csPlot==34)=7;
    csPlot(csPlot==36)=8;
    csSubPlot(csPlot==38)=9;

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
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,1200]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(3,1,1);
    hold on;
    
    sub1=surf(plotTime,plotASL./1000,plotDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    ylim([-0.2,ylimUpper]);
    ylabel('Altitude (km)');
    xlim([plotTime(1),plotTime(end)]);
    title('Reflectivity')
    grid on

    ax2=subplot(3,1,2);
    surf(plotTime,plotASL./1000,csPlot,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([plotTime(1),plotTime(end)]);
    ax2.Colormap=colmapSC;
    caxis([0.5 9.5]);
    cb=colorbar;
    cb.Ticks=1:9;
    cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
        'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
    grid on
    box on
    title('Echo type')
    
    ax3=subplot(3,1,3);
    
    colMapIn=jet(cloudCount-1);
    % Make order random
    indsCol=randperm(size(colMapIn,1));
    colMapInds=cat(2,indsCol',colMapIn);
    colMapInds=sortrows(colMapInds);
    colMap=cat(1,[0 0 0],colMapInds(:,2:end));
    
    hold on;
    sub2=surf(plotTime,plotASL./1000,plotPuzzle,'edgecolor','none');
    view(2);
    ax3.Colormap=colMap;
    caxis([-0.5 cloudCount-1+0.5])
    ylim([-0.2,ylimUpper]);
    ylabel('Altitude (km)');
    xlim([plotTime(1),plotTime(end)]);
    title('Cloud Puzzle')
    grid on
    colorbar
    
    formatOut = 'yyyymmdd_HHMM';    
    set(gcf,'PaperPositionMode','auto')
    print([figdir,project,'_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_cloudPuzzle.png'],'-dpng','-r0');
    
end