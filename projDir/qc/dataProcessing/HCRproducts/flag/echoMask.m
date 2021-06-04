function [echoMask antStat] = echoMask(data)

%% Antenna status (down=0, up=1, pointing=2, scanning=3, transition=4)

antStat=nan(size(data.time));

downInd=find(movmean(data.elevation,10)<-88);
antStat(downInd)=0;
upInd=find(movmean(data.elevation,10)>88);
antStat(upInd)=1;

downInd=find(movmean(data.elevation,500)<-88 & isnan(antStat));
antStat(downInd)=0;
upInd=find(movmean(data.elevation,500)>88 & isnan(antStat));
antStat(upInd)=1;

antDiff=diff(data.elevation);

scan=zeros(size(data.time));
scanInd=find(movmean(abs(antDiff),60)>0.1);
scan(scanInd)=1;

transInd=find(abs(antDiff)>2);
antStat(transInd)=4;

scan(transInd)=0;
ssdiff=diff(scan);
ones1=find(ssdiff==1);
minones1=find(ssdiff==-1);

if ~isempty(ones1) | ~isempty(minones1)
    if ones1(1)>minones1(1)
        ones1=cat(2,1,ones1);
    end
    if length(ones1)~=length(minones1)
        minones1=cat(2,minones1,length(scan));
    end
    
    for ii=1:length(ones1)
        calLen=minones1(ii)-ones1(ii);
        if calLen<60
            scan(ones1(ii)+1:minones1(ii))=0;
        end
    end
    
    antStat(scan==1)=3;
end

antStat(isnan(antStat))=2; % pointing

echoMask=nan(size(data.DBZ));

%% Find bang (6)
echoMask(1:17,:)=6;

%% Find water surface (7),land surface (8), and below surface (9)
[linInd rowInd rangeToSurf] = hcrSurfInds(data);

dbzRatio=nan(size(rangeToSurf));
linAboveSurf=nan(size(rangeToSurf));
linInSurfDBZ=nan(size(rangeToSurf));

rowInd2=rowInd;

for ii=1:size(rowInd,2)
    if ~isnan(rowInd(ii)) & echoMask(rowInd(ii),ii)~=2
        surfUp=data.DBZ(rowInd(ii)-7:rowInd(ii)-1,ii);
        if abs(max(surfUp)-data.DBZ(rowInd(ii),ii))<1
            rowInd2(ii)=rowInd(ii)-1;
        end
        surfStart=5;
        if isfield(data,'TOPO')
            topoRay=data.TOPO(ii);
            if topoRay~=0
                surfStart=10;
            end
        end
        
        if ~isfield(data,'TOPO') | data.TOPO(ii)==0 % Water
            echoMask(rowInd2(ii)-surfStart:min(rowInd2(ii)+7,size(echoMask,1)),ii)=7;
        else % Land
            echoMask(rowInd2(ii)-surfStart:min(rowInd2(ii)+7,size(echoMask,1)),ii)=8;
        end
        % Below surface
        echoMask(min(rowInd2(ii)+8,size(echoMask,1)):end,ii)=9;
        
        % Find ratio of surface to above surface data
        aboveSurf=data.DBZ(18:rowInd2(ii)-surfStart-1,ii);
        linAboveSurf(ii)=sum(10.^(aboveSurf./10),'omitnan');
        inSurfDBZ=data.DBZ(rowInd2(ii)-surfStart:min(rowInd2(ii)+7,size(echoMask,1)),ii);
        linInSurfDBZ(ii)=sum(10.^(inSurfDBZ./10),'omitnan');
        
        dbzRatio(ii)=linInSurfDBZ(ii)/linAboveSurf(ii);
     
    end
end

%% Out of range (5)
if isfield(data,'TOPO')
    compAlt=data.altitude-data.TOPO;
else
    compAlt=data.altitude;
end

theoSurf=compAlt./cosd(data.elevation+90);
theoSurf(data.elevation>0)=nan;

surfDist=theoSurf-data.range(end,:);
surfDist(surfDist<0)=nan;
rangeSpace=data.range(2,1)-data.range(1,1);

outRangePix=ceil(surfDist./rangeSpace);

for ii=1:size(echoMask,2)
    if ~isnan(outRangePix(ii))
        echoMask(18:min(outRangePix(ii),size(echoMask,1)),ii)=5;
    end
end

%% Backlobe echo (4)

% Initiate mask
blMask=zeros(size(data.range)); % When using thresholds

% DBZ, SW threshold
if isfield(data,'WIDTH_CORR')
    data.WIDTH=data.WIDTH_CORR;
end

blMask(data.DBZ<-20 & data.WIDTH>1.4)=1;

% Remove small areas
blMask=bwareaopen(blMask,10);

% Fill holes
blMask=imfill(blMask,'holes');

% Only within right altitude
if isfield(data,'topo')
    rightAlt=(data.altitude-data.topo)*2;
elseif isfield(data,'TOPO')
    rightAlt=(data.altitude-data.TOPO)*2;
else
    rightAlt=data.altitude*2;
end
altMat=repmat(rightAlt,size(data.range,1),1);
% Lower limit
blMask(data.range<(altMat-600))=0;
% Upper limit
blMask(data.range>(altMat+600))=0;

% Only when scanning up
blMask(:,find(data.elevation<0))=0;

% Not in bang
blMask(echoMask==6)=0;

echoMask(blMask==1)=4;

%% Transition (11)

echoMask(:,find(antStat==4))=11;

%% Missing data (12)
bangData=data.DBZ(1:14,:);
bangData(bangData<-10)=nan;
bangMask=zeros(size(bangData));
bangMask(isnan(bangData))=1;
sumNan=sum(bangMask,1,'omitnan');

missInds=find(sumNan>13);
echoMask(:,missInds)=12;

%% NS cal (10)
firstGate=data.DBMVC(1,:);
nscalInds=find(firstGate>-96 & firstGate<-87);
nscalMask=zeros(size(firstGate));
nscalMask(nscalInds)=1;
nscalMask(~isnan(outRangePix))=0;

echoMask(:,find(nscalMask==1))=10;

%% Cloud echo (1)

refl=data.DBZ;
refl(echoMask>0)=nan;

goodInds=find(~isnan(refl));
echoMask(goodInds)=1;

%% Extinct (3)
maskTemp=flipud(echoMask);

reflTemp=flipud(data.DBZ);
reflTemp(maskTemp>1)=nan;

noSurfInds=zeros(size(dbzRatio));
noSurfInds(dbzRatio<0.5)=1;
noSurfInds(isnan(dbzRatio))=1;

% Only when surface is ocean
oceanInds=any(echoMask==7,1);
oceanLowInds=find(linInSurfDBZ<10000 & oceanInds==1);
noSurfInds(oceanLowInds)=1;

% Not when last range is higher than topo
lastAsl=data.asl(end,:);
diffLast=lastAsl-data.TOPO;
noSurfInds(lastAsl-data.TOPO>0)=0;

% Not when last index is cloud
lowest=echoMask(end,:);
noSurfInds(lowest==1)=0;

noSurfInds(linAboveSurf<300)=0;
noSurfInds(~isnan(outRangePix))=0;
noSurfInds(data.elevation>0)=0;
noSurfCols=find(noSurfInds==1);

for ii=1:length(noSurfCols)
    reflCol=reflTemp(:,noSurfCols(ii));
    nanCol=find(isnan(reflCol));
    diffCol=diff(nanCol);
    bigInd=min(find(diffCol>10));
    firstInd=nanCol(bigInd);
    maskTemp(1:firstInd,noSurfCols(ii))=3;
end
echoMask=flipud(maskTemp);

%% Speckle
specCut=100;

reflTemp=data.DBZ;
reflTemp(echoMask>1)=nan;

maskTemp=zeros(size(reflTemp));
maskTemp(~isnan(reflTemp))=1;

CC = bwconncomp(maskTemp);

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=specCut
         echoMask(area)=2;
    end
end

end
