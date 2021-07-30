% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/home/romatsch/plots/HCR/meltingLayer/paper/'];

indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/otrec/10hz/';

% Loop through cases
startTime=datetime(2019,9,21,12,26,0);
endTime=datetime(2019,9,21,12,31,0);

%% Load data

disp('Loading data ...');

data=[];

data.DBZ=[];
data.LDR=[];
data.VEL_CORR=[];
data.TEMP=[];
data.WIDTH=[];
data.FLAG=[];
data.TOPO=[];

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

%% Find melting layer
data.dbzMasked=data.DBZ;
data.dbzMasked(data.FLAG>1)=nan;

[meltLayer iceLayer offset]=f_meltLayer(data,-170);
zeroAdjustMeters=-170;

debugFig=1;

%% Find zero degree altitude

disp('Searching zero degree altitude ...');

oneGate=data.range(2)-data.range(1);
zeroAdjustGates=round(zeroAdjustMeters/oneGate);

tempTemp=data.TEMP;
%tempTemp(1:16,:)=nan;
oddAngles=find(data.elevation>-70 & data.elevation<70);
tempTemp(:,oddAngles)=nan;
signTemp=sign(tempTemp);
zeroDeg=diff(signTemp,1);

clear tempTemp signTemp

zeroDeg(isnan(zeroDeg))=0;
zeroDeg=cat(1,zeroDeg,zeros(size(data.time)));
zeroDeg(zeroDeg~=0)=1;

zeroSum=sum(zeroDeg,1);
tempDataY=find(zeroSum~=0);
zeroSum=zeroSum(tempDataY);
zeroDeg=zeroDeg(:,tempDataY);

zeroAlt=nan(size(zeroDeg));
aslTemp=data.asl(:,tempDataY);
zeroAlt(zeroDeg==1)=aslTemp(zeroDeg==1);

%% Connect zero degree layers

disp('Connecting zero degree layers ...');

% Find how many melting layers there are and connect the right ones
zeroTemp=zeroDeg;

numZero=sum(zeroTemp,1,'omitnan');

clear zeroTemp

% Connect layers
layerInds=nan(1,length(tempDataY));
layerAlts=nan(1,length(tempDataY));

for ii=1:length(tempDataY)
    if numZero==0
        continue
    end
    colInds=find(zeroDeg(:,ii)==1);
    colAltsAll=zeroAlt(:,ii);
    colAlts=colAltsAll(colInds);
    
    if ii>1
        prevAlts=layerAlts(:,ii-1);
    end
    
    % Go through all layers
    for jj=1:numZero(ii)
        % First data points or after nans
        if ii==1 | numZero(ii-1)==0
            if ii==1 & jj==1
                layerInds(jj,ii)=colInds(jj);
                layerAlts(jj,ii)=colAlts(jj);
            else
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        else % Find closest altitude
            zeroAltDiffs=abs(prevAlts-colAlts(jj));
            zeroAltDiffs(isnan(zeroAltDiffs))=100000;
            minDiff=min(zeroAltDiffs);
            if minDiff<50 | (numZero(ii)==1 & numZero(ii-1)==1)
                diffInd=find(zeroAltDiffs==minDiff);
                layerInds(diffInd,ii)=colInds(jj);
                layerAlts(diffInd,ii)=colAlts(jj);
            elseif minDiff~=100000;
                % Add new row
                layerInds=cat(1,layerInds,nan(size(tempDataY)));
                layerAlts=cat(1,layerAlts,nan(size(tempDataY)));
                layerInds(size(layerInds,1),ii)=colInds(jj);
                layerAlts(size(layerInds,1),ii)=colAlts(jj);
            end
        end
    end
end

clear zeroAlt zeroDeg

%% Add adjusted layers

disp('Adjusting zero degree layers ...');
layerIndsAdj=nan(size(layerInds));
layerAltsAdj=nan(size(layerAlts));

adjustMeters=zeroAdjustGates*oneGate;
elevTemp=data.elevation(tempDataY);

for kk=1:length(tempDataY)
    if elevTemp(kk)<0
        layerIndsAdj(:,kk)=layerInds(:,kk)-zeroAdjustGates;
    else
        layerIndsAdj(:,kk)=layerInds(:,kk)+zeroAdjustGates;
    end
    layerAltsAdj(:,kk)=layerAlts(:,kk)+adjustMeters;
end

%% Remove data that is not suitable
LDRdata=data.LDR(:,tempDataY);
tempFlag=data.FLAG(:,tempDataY);
%data=rmfield(data,'FLAG');
LDRdata(tempFlag>1)=nan;
LDRdata(LDRdata<-16 | LDRdata>-7)=nan;
LDRdata(1:20,:)=nan;

VELdata=data.VEL_CORR(:,tempDataY);
VELdata(tempFlag>1)=nan;
VELdata(1:20,:)=nan;

clear tempFlag

%% Tightened backlobe
% Initiate mask
blMask=zeros(size(aslTemp));

tempDBZ=data.DBZ(:,tempDataY);
tempWIDTH=data.WIDTH(:,tempDataY);
blMask(tempDBZ<-18 & tempWIDTH>1)=1;

clear tempWIDTH tempDBZ

% Remove small areas
%blMask(:,find(elevTemp<0))=1;
blMask=bwareaopen(blMask,10);

% Fill holes
blMask=imfill(blMask,'holes');

% Only when scanning up
blMask(:,find(elevTemp<0))=0;

% Only within right altitude
rightAlt=(data.altitude(tempDataY)-data.TOPO(tempDataY))*2;

altMat=repmat(rightAlt,size(data.asl(:,tempDataY),1),1);
% Lower limit
blMask(data.asl(:,tempDataY)<(altMat-600))=0;
% Upper limit
blMask(data.asl(:,tempDataY)>(altMat+600))=0;

LDRdata(blMask==1)=nan;
VELdata(blMask==1)=nan;

clear blMask

%% Find altitude of bright band

disp('Searching for melting layer ...');

BBaltAll={};
timeIall={};
BBaltInterpAll={};
BBaltZeroAll={};
timeIZeroAll={};

% Transition length should be 600 seconds
resol=seconds(median(diff(data.time)));
transLength=1/resol*600;

BBaltRawAll=cell(1,size(layerAltsAdj,1));

kk=1;
timeInds=find(~isnan(layerAltsAdj(kk,:)) & layerIndsAdj(kk,:)>0);
layerAltsTemp=layerAltsAdj(kk,timeInds);
rowInds=layerIndsAdj(kk,timeInds);

% LDR
% Find maximum LDR level in area around zero degree altitude
maxLevelLDR=nan(size(rowInds));
for ii=1:length(rowInds)
    vertColLDR=LDRdata(:,timeInds(ii));
    vertColLDR(rowInds(ii)+18:end)=nan;
    vertColLDR(1:max([rowInds(ii)-18,1]))=nan;
    
    % Check if all nan
    if min(isnan(vertColLDR))==0
        %vertColData=vertColLDR(~isnan(vertColLDR));
        maxLevelLDR(ii)=min(find(vertColLDR==nanmax(vertColLDR)));
    end
end

ASLlayer=aslTemp(:,timeInds);

if min(isnan(maxLevelLDR))==0
    
    % LDR
    colIndsLDR=1:1:size(ASLlayer,2);
    linearIndLDR = sub2ind(size(ASLlayer), maxLevelLDR, colIndsLDR);
    linearIndLDR(isnan(maxLevelLDR))=[];
    
    % Raw altitude of melting layer
    LDRaltRawIn=ASLlayer(linearIndLDR);
    LDRaltRaw=nan(size(rowInds));
    LDRaltRaw(find(~isnan(maxLevelLDR)))=LDRaltRawIn;
    % Mean altitude of meltig layer
    LDRalt=movmedian(LDRaltRaw,300,'omitnan');
    LDRalt(isnan(LDRaltRaw))=nan;
    % Standard deviation
    LDRaltS=movstd(LDRaltRaw,300,'omitnan');
    LDRaltS(isnan(LDRaltRaw))=nan;
    % Remove data with too much std
    LDRaltRaw(LDRaltS>100)=nan;
    % Distance between raw and mean altitude
    LDRloc=abs(LDRaltRaw-LDRalt);
    
    % Remove data where distance is more than 50 m
    LDRaltRaw(LDRloc>50)=nan;
    
    % Adjust zero degree layer
    zeroDist=LDRaltRaw-layerAltsTemp;
    
    % Remove small data stretches
    zeroMask=zeros(size(zeroDist));
    zeroMask(~isnan(zeroDist))=1;
    
    zeroMask=bwareaopen(zeroMask,5);
    zeroDist(zeroMask==0)=nan;
    
    zeroDistSmooth=movmean(zeroDist,3000,'omitnan');
    zeroDistSmooth(isnan(zeroDist))=nan;
    
    noNanDist=fillmissing(zeroDistSmooth,'linear','EndValues','nearest');
    if sum(~isnan(noNanDist))>0
        layerAltsTemp=layerAltsTemp+noNanDist;
        
        distInds=round(noNanDist./oneGate);
        rowInds(elevTemp(timeInds)>0)=rowInds(elevTemp(timeInds)>0)-distInds(elevTemp(timeInds)>0);
        rowInds(elevTemp(timeInds)<=0)=rowInds(elevTemp(timeInds)<=0)+distInds(elevTemp(timeInds)<=0);
    end
end

% VEL

% Find vel melting layer level
maxLevelVEL=nan(size(rowInds));
velMasked=VELdata;
velMasked(:,elevTemp<0)=-velMasked(:,elevTemp<0);

velTime=velMasked(:,timeInds);
for ii=1:length(rowInds)
    velTime(rowInds(ii)+50:end,ii)=nan;
    velTime(1:max([rowInds(ii)-50,1]),ii)=nan;
end

% Find discontinuity in vel field
% Smooth in range direction
velSmooth=movmedian(velTime,10,1);

% Calculate velocity steps
[velSteps,S1,S2] = ischange(velSmooth,1,'MaxNumChanges',1);

% Find altitude of step
stepInLin=find(velSteps==1);
[stepInR,stepInC]=ind2sub(size(velSmooth),stepInLin);

maxLevelVEL=nan(1,length(timeInds));
maxLevelVEL(stepInC)=stepInR;

% Remove data that doesn't cut it

% VEL
colIndsVEL=1:1:size(ASLlayer,2);
linearIndVEL = sub2ind(size(ASLlayer), maxLevelVEL, colIndsVEL);
linearIndVEL(isnan(maxLevelVEL))=[];

% Raw altitude of melting layer
VELaltRawIn=ASLlayer(linearIndVEL);
VELaltRaw=nan(size(rowInds));
VELaltRaw(find(~isnan(maxLevelVEL)))=VELaltRawIn;

% Difference between steps
diffS1=nan(1,size(velSteps,2));
udS1=nan(2,size(velSteps,2));

for ii=1:size(velSteps,2)
    uS1=unique(S1(:,ii),'stable');
    uS1(isnan(uS1))=[];
    
    if elevTemp(ii)>0
        uS1=flipud(uS1);
    end
    
    if length(uS1)==2
        diffS1(ii)=uS1(2)-uS1(1);
        udS1(:,ii)=uS1;
    end
    
end

VELaltRawOrig=VELaltRaw;

VELaltRaw(movmean(udS1(1,:),50,'omitnan')>0)=nan;
        VELaltRaw(movmean(udS1(2,:),50,'omitnan')>0)=nan;
        
updraft=cat(2,find(movmean(udS1(1,:),50,'omitnan')>0),...
    find(movmean(udS1(2,:),50,'omitnan')>0));
updraft=unique(updraft);

VELaltRaw(movmean(diffS1,50,'omitnan')>-0.7)=nan;
diffAB=find(movmean(diffS1,50,'omitnan')>-0.7);

% Compare with zero degree layer
VELzeroDiff=VELaltRaw-layerAltsTemp;
VELaltRaw(abs(VELzeroDiff)>200)=nan;
diffLDR=find(abs(VELzeroDiff)>200);

% Mean altitude of meltig layer
VELalt=movmedian(VELaltRaw,1000,'omitnan');
VELalt(isnan(VELaltRaw))=nan;

% Distance between raw and mean altitude
VELloc=abs(VELaltRaw-VELalt);
VELaltRaw(VELloc>100)=nan;
diffMovmean=find(VELloc>100);

% Standard deviation
VELaltS=movstd(VELaltRaw,300,'omitnan');
VELaltS(isnan(VELaltRaw))=nan;
VELaltRaw(VELaltS>35)=nan;
diffSTD=find(VELaltS>35);

%% Plot
close all

%fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1400,800]);
wi=5;
hi=6.5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');
colmap=jet;
colormap(flipud(colmap));
s1=subplot(2,1,1)
hold on
surf(data.time,data.asl./1000,velSmooth,'edgecolor','none');
view(2)
hcb=colorbar('north');
caxis([-5 5])
plot(data.time,VELaltRawOrig./1000,'-k','linewidth',1.5)
xlim([startTime endTime]);
ylim([3.7 5.9])
grid on
title('(a) Smoothed radial velocity (m s^{-1})');
ylabel('Altitude (km)');
s1.Position=[0.12 0.57 0.83 0.38];
hcb.Position=[0.56 0.9 0.35 0.025];

s2=subplot(2,1,2);
hold on
plot(data.time,VELaltRawOrig./1000,'-k','linewidth',1.5)
scatter(data.time(updraft),VELaltRawOrig(updraft)./1000,'filled');
scatter(data.time(diffLDR),VELaltRawOrig(diffLDR)./1000,'filled');
scatter(data.time(diffMovmean),VELaltRawOrig(diffMovmean)./1000,'filled');
scatter(data.time(diffSTD),VELaltRawOrig(diffSTD)./1000,'filled');
xlim([startTime endTime]);
ylim([3.7 5.9])
grid on
legend('Discontinuity','Upward motion','Distance to LDR melt. alt.',...
    'Distance to moving avg.','Standard deviation');
title('(b) Discontinuity and eliminated data');
ylabel('Altitude (km)');
s2.Position=[0.12 0.08 0.83 0.38];

print([figdir,'velCleanup'],'-dpng','-r0');
