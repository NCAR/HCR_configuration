% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2019,8,11,17,40,0);
endTime=datetime(2019,8,11,17,55,0);

% startTime=datetime(2019,10,2,15,0,0);
% endTime=datetime(2019,10,2,15,59,0);

plotTest=0;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

startThresh=-20; % starting threshold for cloud segmentation
ylimits=[0 15];

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudPuzzle/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

indir=HCRdir(project,quality,freqData);

%% Load data

disp('Loading data ...');

data.DBZ=[];
data.TEMP=[];
data.FLAG=[];

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

cloudPuzzleOut=nan(size(data.DBZ));

refl=data.DBZ;
refl(data.FLAG>1)=nan;

%% Handle missing and NS cal

disp('Filling missing and NS cal ...');

refl = fillMissingNScal(refl,data);

%% Identify contiguous clouds that are too small

disp('Identifying small clouds ...');

%reflLarge=refl;

reflMask=zeros(size(refl));
reflMask(~isnan(refl))=1;

pixCut=5000;
CC = bwconncomp(reflMask);

cloudNumOrig=nan(size(data.DBZ));
countCloud=1;

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=pixCut
        %reflLarge(area)=nan;
        cloudPuzzleOut(area)=0;
    else
        cloudNumOrig(area)=countCloud;
        countCloud=countCloud+1;
    end
end

%% Add extinct back in

disp('Filling extinct echo ...');

[cloudNum reflExt]=fillExtinct(data,cloudNumOrig,refl);

%% Smooth with convolution

% Create mask for convolution
radius=5;
numPix=radius*2+1;
[rr cc] = meshgrid(1:numPix);
cirMask = sqrt((rr-(radius+1)).^2+(cc-(radius+1)).^2)<=radius;
cirMask=double(cirMask);

% Normalize
cirMask=cirMask./(sum(reshape(cirMask,1,[])));

% Convolution
reflConv=nanconv(reflExt,cirMask);

%% Split up individual clouds

disp('Splitting clouds ...');

numMax=max(reshape(cloudNum,1,[]),[],'omitnan');

thresh12=[-25 0];

cloudCount=1;

for ii=1:numMax
    
    cloudInds=find(cloudNum==ii);
    
    if length(cloudInds)>100000
             
        reflMapBig=nan(size(cloudNum));
        reflMapBig(cloudInds)=reflConv(cloudInds);
        
        [clR clC]=ind2sub(size(cloudNum),cloudInds);
        
        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
        aslMap=data.asl(min(clR):max(clR),min(clC):max(clC));
        timeMap=data.time(min(clC):max(clC));
        
        % Zero padding
        reflPadded=cat(1,nan(10,size(reflMap,2)),reflMap,nan(10,size(reflMap,2)));
        reflPadded=cat(2,nan(size(reflPadded,1),10),reflPadded,nan(size(reflPadded,1),10));
                     
        % BW mask with filled holes
        BW=zeros(size(reflPadded));
        BW(~isnan(reflPadded))=1;
        
        BW2 = imfill(BW,'holes');
                 
        % Distance
        D = -bwdist(~BW2);
               
        % Thresholding
        maskBig = thresholdMask(reflPadded,startThresh,500);
                
        newMin = imimposemin(D,maskBig);
        waterShed = watershed(newMin);
        
        waterMasked=joinCloudParts(waterShed,BW2);
        
        maskJoined=zeros(size(BW2));
        maskJoined(waterMasked>0)=1;
        
        maskJoined=maskJoined(11:end-10,11:end-10);
        
        maskBack=zeros(size(reflMapBig));
        maskBack(min(clR):max(clR),min(clC):max(clC))=maskJoined;
                
        uniqueClouds=bwconncomp(maskBack);
        
        for kk=1:uniqueClouds.NumObjects
            areaU=uniqueClouds.PixelIdxList{kk};
            cloudPuzzleOut(areaU)=cloudCount;
            cloudCount=cloudCount+1;
        end
        
        % Plot sub plot with individual cloud
        if plotTest
            close all
            cloudOut=nan(size(maskJoined));
            
            uniqueClouds2=bwconncomp(maskJoined);
            
            for kk=1:uniqueClouds2.NumObjects
                areaU=uniqueClouds2.PixelIdxList{kk};
                cloudOut(areaU)=kk;
            end
            
            timeMat=repmat(timeMap,size(reflMap,1),1);
            
            fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
            
            s1=subplot(2,1,1);
            hold on
            sub1=surf(timeMap,aslMap./1000,reflMap,'edgecolor','none');
            view(2);
            sub1=colMapDBZ(sub1);
            ylabel('Altitude (km)');
            xlim([timeMap(1),timeMap(end)]);
            title('Reflectivity')
            grid on
            pos1=s1.Position;
            
            s2=subplot(2,1,2);
            hold on
            sub2=surf(timeMap,aslMap./1000,cloudOut,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            xlim([timeMap(1),timeMap(end)]);
            grid on
            pos2=s2.Position;
            s2.Position=[pos1(1),pos2(2),pos1(3),pos1(4)];
        end
    else
        cloudPuzzleOut(cloudInds)=cloudCount;
        cloudCount=cloudCount+1;
    end
end
%% Plot

disp('Plotting ...');

close all

fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);

%%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1=subplot(2,1,1);
hold on;

sub1=surf(data.time,data.asl./1000,refl,'edgecolor','none');
view(2);
sub1=colMapDBZ(sub1);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity')
grid on

ax2=subplot(2,1,2);

colMap=lines;
colMap=cat(1,[0 0 0],colMap);

hold on;
sub2=surf(data.time,data.asl./1000,cloudPuzzleOut,'edgecolor','none');
view(2);
ax2.Colormap=colMap;
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Cloud Puzzle')
grid on

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_cloudPuzzle.png'],'-dpng','-r0');

