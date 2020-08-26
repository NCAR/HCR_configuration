% Cloud classification algorithm
% Small clouds or not classified because e.g. plane in cloud (0)
% Low clouds: Deep (1), Ns (2), Cu (3), Sc (4), St (5)
% Middle clouds: As (6), Ac (7)
% High clouds (8)

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,10,2,15,0,0);
endTime=datetime(2019,10,2,15,59,0);

%getFreezeL=0;

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
data.TOPO=[];
data.TEMP=[];

% if getFreezeL
%     data.LDR=[];
%     data.VEL_CORR=[];
%     data.WIDTH=[];
% else
%     data.FREEZING_LEVEL=[];
% end

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
% 
% if getFreezeL
%     data.FREEZING_LEVEL=f_meltLayer(data,200);
% end

%% Cloud puzzle

data.cloudPuzzle=f_cloudPuzzle_radial(data);

%% Surface refl
[linInd rowInd rangeToSurf] = hcrSurfInds(data);
data.surfRefl=data.DBZ(linInd);

%% Loop through clouds

puzzleReplace=data.cloudPuzzle;
puzzleReplace(isnan(puzzleReplace))=-99;

cloudNums=unique(puzzleReplace);
cloudNums(cloudNums==0 | cloudNums==-99)=[];

% Output
cloudClass=nan(size(data.DBZ));

for ii=1:length(cloudNums)
    [rowInd colInd]=find(data.cloudPuzzle==cloudNums(ii));
    wholeInd=find(data.cloudPuzzle==cloudNums(ii));
    
    disp(['Identifying cloud ',num2str(ii),' of ',num2str(length(cloudNums)),' ...']);
    
    %% Cut out columns with cloud data
    allVars=fields(data);
    
    for jj=1:size(allVars,1)
        varIn=data.(allVars{jj});
        startVar=min(colInd);
        if startVar~=1
            startVar=startVar-1;
        end
        endVar=max(colInd);
        if endVar~=length(data.time)
            endVar=endVar+1;
        end
        if min(size(varIn))==1
            dataCut.(allVars{jj})=varIn(startVar:endVar);
        else
            dataCut.(allVars{jj})=varIn(:,startVar:endVar);
        end
    end
    
    dataCut.puzzleOne=dataCut.cloudPuzzle;
    dataCut.puzzleOne(dataCut.puzzleOne~=cloudNums(ii))=nan;
        
    dataCut.asl=dataCut.asl./1000; % In km
    dataCut.asl(dataCut.puzzleOne~=cloudNums(ii))=nan;
    dataCut.TEMP(dataCut.puzzleOne~=cloudNums(ii))=nan;
    dataCut.DBZ(dataCut.puzzleOne~=cloudNums(ii))=nan;
    dataCut.FLAG(dataCut.FLAG==8)=7;
    dataCut.flagCensored=dataCut.FLAG;
    dataCut.flagCensored(dataCut.puzzleOne~=cloudNums(ii))=nan;
    
    %% Calculate cloud parameters
    cloudParams=calcCloudParams(dataCut);
    
    %% Cloud classifier
    % Remove cloudParams if they are nan so they cannot be used. Then
    % try statement below will go into catch.
    fieldParams=fields(cloudParams);
    
    for ii=1:length(fieldParams)
        if isnan(cloudParams.(fieldParams{ii})) & ~strcmp(fieldParams{ii},'max10dbzAgl')
            cloudParams=rmfield(cloudParams,(fieldParams{ii}));
        end
    end
    
    % Cloud classification
    cloudFlag=[];
    try
        if cloudParams.meanMinAgl>10 & ~cloudParams.precip % High cloud (8)
            cloudFlag=8;
        elseif (cloudParams.meanMaxAgl>2.5 & cloudParams.precip) % Precipitating clouds: Deep (1), Ns (2), Cu (3), Sc (4), St (5), Ac (7)
            cloudFlag=precipCloudClass(cloudParams);
        elseif (cloudParams.meanMaxReflTemp<-23 & cloudParams.meanMaxRefl<-3 & ...
                cloudParams.meanMaxReflAgl>5 & cloudParams.meanMinAgl>5)
            cloudFlag=highCloudClass(cloudParams); % High clouds: Deep (1), Cu (3), As (6), High (8)
        elseif (cloudParams.meanMaxReflTemp>-15 & cloudParams.meanMaxReflAgl<2) ...
                | cloudParams.meanMinAgl<1.5 % Low clouds: Deep (1), Ns (2), Cu (3), Sc (4), St (5), As (6), Ac (7)
            cloudFlag=lowCloudClass(cloudParams);
        else
            cloudFlag=middleCloudClass(cloudParams); % Middle clouds: Ns (2), Cu (3), Sc (4), St (5), As (6), Ac (7)
        end
        
        cloudClass(wholeInd)=cloudFlag;
    catch
        cloudClass(wholeInd)=0;
    end
end

cloudClass(data.cloudPuzzle==0)=0; % Small, unclassified echos

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
ax2.Colormap=lines;
hold on;
sub3=surf(data.time,data.asl./1000,data.cloudPuzzle,'edgecolor','none');
view(2);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Current cloud')
grid on

ax3=subplot(3,1,3);
colmap=jet(8);
ax3.Colormap=cat(1,[0 0 0],colmap);
hold on;
sub3=surf(data.time,data.asl./1000,cloudClass,'edgecolor','none');
view(2);
caxis([0 9]);
colorbar('Ticks',(0.5:1:9.5),'YTickLabel',{'N/A','Deep','Ns','Cu','Sc','St','As','Ac','High'});
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Cloud classification')
grid on

%     formatOut = 'yyyymmdd_HHMM';
%     set(gcf,'PaperPositionMode','auto')
%     print([figdir,'cloudID',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
