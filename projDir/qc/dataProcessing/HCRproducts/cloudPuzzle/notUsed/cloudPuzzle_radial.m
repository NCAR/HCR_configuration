% Analyze HCR clouds

clear all;
close all;

plotTest=1;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

if strcmp(project,'otrec')
    ylimits=[-0.2 15];
else
    ylimits=[-0.2 10];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudPuzzle/'];
figdir=['/home/romatsch/plots/HCR/cloudPuzzle/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

% Loop through cases
casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/cloudPuzzle_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=2:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
    data.DBZ=[];
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
    
    reflMask=zeros(size(refl));
    reflMask(~isnan(refl))=1;
    
    pixCut=5000;
    CC = bwconncomp(reflMask);
    
    cloudNumOrig=nan(size(data.DBZ));
    countCloud=1;
    
    for ii=1:CC.NumObjects
        area=CC.PixelIdxList{ii};
        if length(area)<=pixCut
            cloudPuzzleOut(area)=0;
        else
            cloudNumOrig(area)=countCloud;
            countCloud=countCloud+1;
        end
    end
    
    %% Add extinct back in
    
    disp('Filling extinct echo ...');
    
    [cloudNum,reflExt]=fillExtinct(data,cloudNumOrig,refl);
        
    %% Split up individual clouds
    
    disp('Splitting clouds ...');
    
    numMax=max(reshape(cloudNum,1,[]),[],'omitnan');
    
    cloudCount=1;
    
    for ii=1:numMax
        
        cloudInds=find(cloudNum==ii);
        
        if length(cloudInds)>100000
            
            reflMapBig=nan(size(cloudNum));
            reflMapBig(cloudInds)=reflExt(cloudInds);
            
            [clR,clC]=ind2sub(size(cloudNum),cloudInds);
            
            reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
            
            % Zero padding
            reflPadded=cat(1,nan(10,size(reflMap,2)),reflMap,nan(10,size(reflMap,2)));
            reflPadded=cat(2,nan(size(reflPadded,1),10),reflPadded,nan(size(reflPadded,1),10));
                         
            BW=zeros(size(reflPadded));
            BW(~isnan(reflPadded))=1;
            
            BW2=imerode(BW, strel('disk', 5));
            BW3 = bwareaopen(BW2,1000);
                        
            % Distance of each cloud pixel from reflectivity threshold mask
            D = bwdist(BW3);
            
            % Watershed is an image segmentation method that looks for
            % ridges and valleys in an image
            waterShed = watershed(D);
            
            % Watershed usually over-segments so we join areas back together
            % that share a large border or where one area is too small
            
            waterCensored=double(waterShed);
            waterCensored(~BW)=nan;
            
            waterMasked=joinCloudParts(waterCensored);
            
            maskJoined=zeros(size(BW2));
            maskJoined(waterMasked>0)=1;
            
            % Reverser zero padding
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
                
                aslMap=data.asl(min(clR):max(clR),min(clC):max(clC));
                timeMap=data.time(min(clC):max(clC));
                
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
    
    cloudPuzzleOut(isnan(reflExt))=nan;
    
    %% Fill in pixels that are not in small areas (i.e. not zero) that have
    % reflectivities but are nan in cloudPuzzle
    
    disp('Filling in final pixels ...');
    allReflMask=zeros(size(reflExt));
    allReflMask(~isnan(reflExt))=1;
    allReflMask(cloudPuzzleOut==0)=0;
    
    puzzleMask=zeros(size(reflExt));
    puzzleMask(cloudPuzzleOut>0)=1;
    
    [oldR oldC]=find(puzzleMask==1);
    [addR addC]=find(puzzleMask==0 & allReflMask==1);
    idx = knnsearch([oldR oldC], [addR addC]);
    nearest_OldValue = cloudPuzzleOut(sub2ind(size(cloudPuzzleOut), oldR(idx), oldC(idx)));
    cloudPuzzleAttached=cloudPuzzleOut;
    cloudPuzzleAttached(sub2ind(size(cloudPuzzleOut), addR, addC))=nearest_OldValue;
    
    % Sometimes areas get attached to wrong area
    cloudPuzzleFinal=cloudPuzzleAttached;
    for jj=1:cloudCount-1
        maskNumber=zeros(size(cloudPuzzleAttached));
        maskNumber(cloudPuzzleAttached==jj)=1;
        individs=bwconncomp(maskNumber);
        if individs.NumObjects>1
            indivClouds=individs.PixelIdxList;
            for ll=1:individs.NumObjects
                if length(indivClouds{ll})<1001
                    maskIndiv=zeros(size(cloudPuzzleAttached));
                    maskIndiv(indivClouds{ll})=1;
                    maskExp=imdilate(maskIndiv, strel('disk', 2));
                    pixExp=cloudPuzzleAttached(find(maskExp==1));
                    pixExp(find(pixExp==0 | pixExp==jj | isnan(pixExp)))=[];
                    if ~isempty(pixExp)
                        pixU=unique(pixExp);
                        cloudPuzzleFinal(indivClouds{ll})=pixU;
                    end
                end
            end
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
    
    colMapIn=jet(cloudCount-1);
    % Make order random
    indsCol=randperm(size(colMapIn,1));
    colMapInds=cat(2,indsCol',colMapIn);
    colMapInds=sortrows(colMapInds);
    colMap=cat(1,[0 0 0],colMapInds(:,2:end));
    
    hold on;
    sub2=surf(data.time,data.asl./1000,cloudPuzzleFinal,'edgecolor','none');
    view(2);
    ax2.Colormap=colMap;
    caxis([-0.5 cloudCount-1+0.5])
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Cloud Puzzle')
    grid on
    colorbar
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,project,'_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_cloudPuzzle.png'],'-dpng','-r0');
    
end