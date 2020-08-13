% Cloud classification algorithm

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,10,2,15,0,0);
endTime=datetime(2019,10,2,15,59,0);

getFreezeL=0;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

ylimits=[-0.2 15];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/cloudClassHCR/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data=[];
data.DBZ=[];
data.FLAG=[];

if getFreezeL
    data.LDR=[];
    data.VEL_CORR=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.TOPO=[];
else
    data.FREEZING_LEVEL=[];
end

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Mask
data.dbzMasked=data.DBZ;
data.dbzMasked(data.FLAG>1)=nan;

%% Get freezing level

if getFreezeL
    data.FREEZING_LEVEL=f_meltLayer(data,200);
end

%% Cloud puzzle

cloudPuzzle=f_cloudPuzzle_radial(data);

%% Loop through clouds

puzzleReplace=cloudPuzzle;
puzzleReplace(isnan(puzzleReplace))=-99;

cloudNums=unique(puzzleReplace);
cloudNums(cloudNums==0 | cloudNums==-99)=[];

for ii=1:length(cloudNums)
    [rowInd colInd]=find(cloudPuzzle==cloudNums(ii));
    %wholeInd=find(cloudPuzzle==cloudNums(ii));
    
    puzzleOne=cloudPuzzle;
    puzzleOne(puzzleOne~=cloudNums(ii))=nan;
    
    puzzleCut=puzzleOne(:,min(colInd):max(colInd));
    
    aslCut=data.asl(:,min(colInd):max(colInd));
    aslCut(isnan(puzzleCut))=nan;
    elevCut=data.elevation(min(colInd):max(colInd));
    altCut=data.altitude(min(colInd):max(colInd));
    
    %% Calculate cloud properties
    
    % min/max asl
    minAsl=nan(1,size(puzzleCut,2));
    maxAsl=nan(1,size(puzzleCut,2));
    
    for jj=1:size(puzzleCut,2)
        aslCol=aslCut(:,jj);
        [minAsl(jj) minIndAsl]=min(aslCol);
        [maxAsl(jj) maxIndAsl]=max(aslCol);
        % Check if flying in cloud
        % Pointing down
        if elevCut(jj)<0 & maxAsl(jj)<10000 & maxIndAsl==18
            maxAsl(jj)=nan;
        end
        % Pointing up
        if elevCut(jj)>=0 & minAsl(jj)<10000 & minIndAsl==18
            minAsl(jj)=nan;
        end
    end
    
    percWanted=0.1;
    
    % Make sure we have enough data and calculate percentiles
    if length(find(~isnan(minAsl)))>length(minAsl)/2
        sortedMin=sort(minAsl,'ascend');
        percIndMin=round(percWanted*length(minAsl));
        minPerc=sortedMin(percIndMin);
    end
    if length(find(~isnan(maxAsl)))>length(maxAsl)/2
        sortedMax=sort(maxAsl,'descend');
        percIndMax=round(percWanted*length(maxAsl));
        maxPerc=sortedMax(percIndMax);
    end
    
    %% Plot
    
    timeMat=repmat(data.time,size(data.DBZ,1),1);
        
    close all
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,900]);
    
    ax1=subplot(3,1,1);
    hold on;
    sub1=surf(data.time,data.asl./1000,data.dbzMasked,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity')
    grid on
    
    %%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax2=subplot(3,1,2);
    hold on;
    sub3=surf(data.time,data.asl./1000,puzzleOne,'edgecolor','none');
    view(2);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Current cloud')
    grid on
    
%     ax3=subplot(3,1,3);
%     ax3.Colormap=jet;
%     hold on;
%     sub3=surf(data.time,data.asl./1000,velMasked,'edgecolor','none');
%     view(2);
%     caxis([-5 5]);
%     colorbar
%     ylim(ylimits);
%     ylabel('Altitude (km)');
%     xlim([data.time(1),data.time(end)]);
%     title('VEL')
%     grid on
    
%     formatOut = 'yyyymmdd_HHMM';
%     set(gcf,'PaperPositionMode','auto')
%     print([figdir,'cloudID',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
    
end