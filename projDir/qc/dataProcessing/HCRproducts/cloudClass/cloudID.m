% Analyze HCR clouds

clear all;
close all;

% startTime=datetime(2018,2,7,18,0,0);
% endTime=datetime(2018,2,8,12,0,0);

startTime=datetime(2019,10,2,15,0,0);
endTime=datetime(2019,10,2,15,59,0);

plotFields=1;
plotWholeFlight=0;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

% Expected bright band altitude. Determines plot zoom.
if strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 15];
elseif strcmp(project,'socrates')
    expBBalt=2;
    ylimits=[-0.2 9];
elseif strcmp(project,'otrec')
    expBBalt=5;
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudID/',project,'/'];

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

cloudIDout=nan(size(data.DBZ));

%% Find zero degree altitude

disp('Searching 0 deg altitude ...');

timeMat=repmat(data.time,size(data.TEMP,1),1);

tempTemp=data.TEMP;
tempTemp(1:16,:)=nan;
signTemp=sign(tempTemp);
zeroDeg=diff(signTemp,1);

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(zeroDeg~=0)=1;

zeroAlt=nan(size(data.asl));
zeroAlt(zeroDeg==1)=data.asl(zeroDeg==1);

zeroInds=find(~isnan(zeroAlt));

%% Mask non cloud

refl=data.DBZ;
refl(data.FLAG>1)=nan;

%% Add extinct back in

flagTemp=data.FLAG;
flagTemp(data.FLAG==8)=7;
surfMask=ones(1,size(flagTemp,2));
surfMask(find(any(data.FLAG==7 | data.FLAG==8,1)))=0;
surfDiff=diff(surfMask);

% Add surface echo back in
holeStart=find(surfDiff==1);
holeStart=holeStart+1;
holeEnd=find(surfDiff==-1);

if holeStart(1)>holeEnd(1)
    holeStart=[1 holeStart];
end
if length(holeStart)~=length(holeEnd)
    holeEnd=[holeEnd,size(flagTemp,2)-1];
end

for ii=1:length(holeStart)
    flagColStart=flagTemp(:,holeStart(ii)-1);
    flagColEnd=flagTemp(:,holeEnd(ii)+1);
    
    minStart=min(find(flagColStart==7));
    minEnd=min(find(flagColEnd==7));
    
    if isempty(minStart)
        minStart=minEnd;
    end
    if isempty(minEnd)
        minEnd=minStart;
    end
    
    newVec=holeStart(ii):holeEnd(ii);
    holeLength=holeEnd(ii)-holeStart(ii);
    
    if holeLength==0
        dataFill=holeStart(ii);
    else        
        x=[holeStart(ii),holeEnd(ii)];
        v=[minStart,minEnd];        
        dataFill=round(interp1(x,v,newVec));
    end
    for jj=1:length(dataFill)
        flagTemp(dataFill(jj):end,newVec(jj))=7;
    end
end

% Fill in extinct echo with median of column above
extInd=find(any(data.FLAG==3,1));

for ii=1:length(extInd)
    reflCol=refl(:,extInd(ii));
    reflMed=median(reflCol,'omitnan');
    
    extCol=flagTemp(:,extInd(ii));
    extData=find(extCol==3);
    refl(extData,extInd(ii))=reflMed;
end

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
reflConv=nanconv(refl,cirMask);

%% Remove contiguous clouds that are too small

reflLarge=reflConv;

reflMask=zeros(size(reflConv));
reflMask(~isnan(reflConv))=1;

pixCut=5000;
CC = bwconncomp(reflMask);

cloudNum=nan(size(data.DBZ));
countCloud=1;

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=pixCut
        reflLarge(area)=nan;
        cloudIDout(area)=1;
    else
        cloudNum(area)=countCloud;
        countCloud=countCloud+1;
    end
end

%% Split up individual clouds
numMax=max(reshape(cloudNum,1,[]),[],'omitnan');

for ii=1:numMax
    
    cloudInds=find(cloudNum==ii);
    
    if length(cloudInds)>100000
        close all;
        cloudRefl=reflLarge(cloudInds);
        
        reflMapBig=nan(size(cloudNum));
        reflMapBig(cloudInds)=cloudRefl;
        
        [clR clC]=ind2sub(size(cloudNum),cloudInds);
        
        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
        aslMap=data.asl(min(clR):max(clR),min(clC):max(clC));
        %zeroAltMap=zeroAlt(min(clR):max(clR),min(clC):max(clC));
        timeMap=data.time(min(clC):max(clC));
        
        %         % Watershed
        bw=zeros(size(reflMap));
        bw(~isnan(reflMap))=1;
        
        %
        %         D = -bwdist(~bw);
        %
        %         mask = imextendedmin(D,50);
        %
        %         D2 = imimposemin(D,mask);
        %         Ld2 = watershed(D2);
        %
        %         I2 = im2double(Ld2);
        %         I2(isnan(reflMap))=nan;
        %
        %         % Find unique cloud parts
        %         noNan=I2;
        %         noNan(isnan(noNan))=-999;
        %         unClouds=unique(noNan);
        %         unClouds(unClouds==-999)=[];
        %         unClouds(unClouds==0)=[];
        %
        %         % Test if unique cloud parts should be separate
        %         if length(unClouds)>1
        %             for aa=1:length(unClouds)
        %                 for bb=aa+1:length(unClouds)
        %                     refl1=reflMap(I2==unClouds(aa));
        %                     refl2=reflMap(I2==unClouds(bb));
        %
        %                     figure
        %                     subplot(2,1,1);
        %                     histogram(refl1);
        %                     subplot(2,1,2);
        %                     histogram(refl2);
        %                     mean(refl1)
        %                     mean(refl2)
        %                 end
        %             end
        %         end
        
        % Number of layers
        
        % Create mask
%         BW2 = imfill(bw,'holes');
%         
        BWext=cat(1,zeros(1,size(bw,2)),bw,zeros(1,size(bw,2)));
        BWdiffOrig=abs(diff(BWext,1,1));
        BWdiffOrig=BWdiffOrig(1:end-1,:);
        BWdiff=BWdiffOrig;
        
        % Connect layer edges

disp('Connecting layer edges ...');

% Find how many melting layers there are and connect the right ones

% Remove contiguous layers that are too small
% if length(timeMap)>9000
%     zeroCut=10;
%     CC = bwconncomp(BWdiff);
%     
%     for ii=1:CC.NumObjects
%         area=CC.PixelIdxList{ii};
%         if length(area)<=zeroCut
%             BWdiff(area)=nan;
%         end
%     end
% end

numZero=sum(BWdiff,1,'omitnan');

% Connect layers
layerInds=nan(1,length(timeMap));
layerAlts=nan(1,length(timeMap));

for aa=1:length(timeMap)
    if numZero(aa)==0
        continue
    end
    colInds=find(BWdiffOrig(:,aa)==1);
    colAltsAll=aslMap(:,aa);
    colAlts=colAltsAll(colInds);
    
    if aa>1
        prevAlts=layerAlts(:,aa-1);
    end
    
    % Go through all layers
    for bb=1:numZero(aa)
        % First data points or after nans
        if aa==1 | numZero(aa-1)==0
            if aa==1 & bb==1
                layerInds(bb,aa)=colInds(bb);
                layerAlts(bb,aa)=colAlts(bb);
            else
                % Add new row
                layerInds=cat(1,layerInds,nan(size(timeMap)));
                layerAlts=cat(1,layerAlts,nan(size(timeMap)));
                layerInds(size(layerInds,1),aa)=colInds(bb);
                layerAlts(size(layerInds,1),aa)=colAlts(bb);
            end
        else % Find closest altitude
            zeroAltDiffs=abs(prevAlts-colAlts(bb));
            zeroAltDiffs(isnan(zeroAltDiffs))=100000;
            minDiff=min(zeroAltDiffs);
            if minDiff<50
                diffInd=find(zeroAltDiffs==minDiff);
                layerInds(diffInd,aa)=colInds(bb);
                layerAlts(diffInd,aa)=colAlts(bb);
            elseif minDiff~=100000;
                % Add new row
                layerInds=cat(1,layerInds,nan(size(timeMap)));
                layerAlts=cat(1,layerAlts,nan(size(timeMap)));
                layerInds(size(layerInds,1),aa)=colInds(bb);
                layerAlts(size(layerInds,1),aa)=colAlts(bb);
            end
        end
    end
end
        
edgeLength=sum(~isnan(layerInds),2);

tooShort=find(edgeLength<10);
layerInds(tooShort,:)=[];
layerAlts(tooShort,:)=[];

        % Metrics per ray
        minAlt=nan(size(timeMap));
        maxAlt=nan(size(timeMap));
        maxRefl=nan(size(timeMap));
        medRefl=nan(size(timeMap));
        numLayers=ones(size(timeMap));
        maxThick=nan(size(timeMap));
        
        for kk=1:size(reflMap,2)
            reflRay=reflMap(:,kk);
            dataRay=find(~isnan(reflRay));
            
            closeGate=aslMap(min(dataRay),kk);
            farGate=aslMap(max(dataRay),kk);
            
            minAlt(kk)=min([closeGate,farGate]);
            maxAlt(kk)=max([closeGate,farGate]);
            
            maxRefl(kk)=max(reflRay,[],'omitnan');
            medRefl(kk)=median(reflRay,'omitnan');
            
%             % Number of layers
%             BWray=BWdiff(:,kk);
%             startBW=find(BWray==1);
%             endBW=find(BWray==-1);
%             
%             thickness=endBW-startBW;
%             maxThick(kk)=max(thickness);
%             largeLayers=length(find(thickness>10));
%             if largeLayers>0
%                 numLayers(kk)=largeLayers;
%             end
        end
               
        maxThickKM=maxThick.*(data.range(2)-data.range(1))./1000;
        
        
        % Plot
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
        
        colIndsAll=1:length(timeMap);
        for aa=1:size(layerAlts,1)
            rowInds=layerInds(aa,:);
            colInds=colIndsAll;
            colInds(isnan(rowInds))=[];
            rowInds(isnan(rowInds))=[];
            
            linPlot=sub2ind(size(aslMap),rowInds,colInds);
            
            rowAlts=layerAlts(aa,:);
            rowAlts(isnan(rowAlts))=[];
            scatter(timeMat(linPlot),rowAlts./1000,10,'filled');
        end
        
        s2=subplot(2,1,2);
        hold on
        plot(timeMap,minAlt./1000,'-g','linewidth',1.5);
        plot(timeMap,maxAlt./1000,'-b','linewidth',1.5);
        %plot(timeMap,numLayers,'-c','linewidth',1.5);
        plot(timeMap,maxThickKM,'-k','linewidth',1.5);
        
        yyaxis right
        plot(timeMap,maxRefl,'-r','linewidth',1.5);
        plot(timeMap,medRefl,'-m','linewidth',1.5);
        xlim([timeMap(1),timeMap(end)]);
        grid on
        pos2=s2.Position;
        s2.Position=[pos1(1),pos2(2),pos1(3),pos1(4)];
        show1=1;
    end
end
%% Plot

close all

if plotFields
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    
    sub1=surf(data.time,data.asl./1000,refl,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(timeMat(zeroInds),data.asl(zeroInds)./1000,10,'k','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and ERA5 freezing level')
    grid on
    
    ax2=subplot(2,1,2);
    hold on;
    sub2=surf(data.time,data.asl./1000,cloudNum,'edgecolor','none');
    view(2);
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    % print([figdir,'refl_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_zeroDegree'],'-dpng','-r0');
end

