function [echoMask,antStat] = echoMask(data)

%% Antenna status (down=1, up=2, pointing=3, scanning=4, transition=5, failure=6)
disp('Working on antenna status ...');

antStat=nan(size(data.time));

% Find transition areas
% Find areas with large standard deviation
movStd=movstd(data.elevation,101);
stdFake=zeros(size(movStd));
stdFake(movStd>5)=movStd(movStd>5);

% Find the peaks with maximum Std
[peakTrans,locs] = findpeaks(stdFake);
% Broaden the peaks
broadPeak=nan(size(stdFake));
broadPeak(locs)=1;
broadPeak=movmean(broadPeak,51,'omitnan');

% Find areas with lots of antenna movement around std peaks
antDiff=diff(data.elevation);
antDiff=cat(2,0,antDiff);
findTrans=abs(movmean(antDiff,11));

% Transision zones
antStat(findTrans>0.2 & broadPeak==1)=5;

% Small scale moving std
movStdSmall=movstd(data.elevation,41);

smoothIndsIn=double(movStdSmall<0.25);
smoothIndsIn(smoothIndsIn==0)=nan;
% Fill holes
smoothIndsLarge=movmean(smoothIndsIn,21,'omitnan');
smoothIndsSmall=movmean(smoothIndsLarge,21,'includenan');
% Remove short
smoothIndsSmall(isnan(smoothIndsSmall))=0;
smoothInds=bwareaopen(smoothIndsSmall,10);

% Remove areas that are classified as transition
smoothInds(antStat==5)=0;

% Loop through smooth areas and classify them as up, down, pointing
smoothAreas=bwconncomp(smoothInds);

% Start loop
for ii=1:smoothAreas.NumObjects
    if median(data.elevation(smoothAreas.PixelIdxList{ii}))>88 % Up pointing
        antStat(smoothAreas.PixelIdxList{ii})=2;
    elseif median(data.elevation(smoothAreas.PixelIdxList{ii}))<-88 % Down pointing
        antStat(smoothAreas.PixelIdxList{ii})=1;
    else
        antStat(smoothAreas.PixelIdxList{ii})=3; % Pointing
    end
end

% Find scanning
scanRate=abs(data.elevation(4:end)-data.elevation(1:end-3));
scanRate=cat(2,0,scanRate,0,0);
scanMask=zeros(size(scanRate));
scanMask(scanRate>0.5)=1;
scanMaskFilt=modefilt(scanMask,[1,99]);

% Again we loop through scanning stretches
scanAreas=bwconncomp(scanMaskFilt);

for ii=1:scanAreas.NumObjects
    if length(scanAreas.PixelIdxList{ii})>500 % If the stretch is too short, set to transition
        antStat(scanAreas.PixelIdxList{ii})=4;
    end
end

% Failure
antStat(isnan(antStat))=6;

% Go through failure stretches. If they are short and next to transition and smooth,
% set to smooth
fmask=antStat==6;

failAreas=bwconncomp(fmask);

% Start loop
for ii=1:failAreas.NumObjects
    if length(failAreas.PixelIdxList{ii})<20
        pixArea=failAreas.PixelIdxList{ii};
        pixBefore=max([1,pixArea(1)-1]);
        pixAfter=min([length(antStat),pixArea(end)+1]);
        antStatBA=[antStat(pixBefore),antStat(pixAfter)];
        if any(antStatBA==5)
            antStatBA(antStatBA==5)=[];
            if antStatBA~=6
                antStat(failAreas.PixelIdxList{ii})=antStatBA;
            end
        end
    end
end

% Go through transition stretches. If they are between failure
% set to failure
tmask=antStat==5;
% Loop through smooth areas and classify them as up, down, pointing
transAreas=bwconncomp(tmask);

% Start loop
for ii=1:transAreas.NumObjects
    pixArea=transAreas.PixelIdxList{ii};
    pixBefore=max([1,pixArea(1)-1]);
    pixAfter=min([length(antStat),pixArea(end)+1]);
    antStatBA=[antStat(pixBefore),antStat(pixAfter)];
    if antStatBA(1)==6 & antStatBA(2)==6
        antStat(transAreas.PixelIdxList{ii})=6;
    end
end

%% Echo mask
disp('Working on flag ...')

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

rightAlt=rightAlt+data.TOPO;
altMat=rightAlt;
%altMat=repmat(rightAlt,size(data.range,1),1);
% Lower limit
blMask(data.asl<(altMat-200))=0; % Default before bug fix was 1000
% Upper limit land
blMask(data.asl>(altMat+600) & data.TOPO>0)=0; % Default before bug fix was 600
% Upper limit ocean
blMask(data.asl>(altMat+200) & data.TOPO==0)=0; % Default before bug fix was 600

% Only when scanning up
blMask(:,find(data.elevation<0))=0;

% Not in bang
blMask(echoMask==6)=0;

echoMask(blMask==1)=4;

%% Missing data (11)
bangData=data.DBZ(1:14,:);
bangData(bangData<-10)=nan;
bangMask=zeros(size(bangData));
bangMask(isnan(bangData))=1;
sumNan=sum(bangMask,1,'omitnan');

missInds=find(sumNan>13);
echoMask(:,missInds)=11;

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

%% Speckle (2)
specCut=100;

reflTemp=data.DBZ;
reflTemp(echoMask>1)=nan;

maskTemp=zeros(size(reflTemp));
maskTemp(~isnan(reflTemp))=1;

CC=bwconncomp(maskTemp);

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=specCut
        echoMask(area)=2;
    end
end

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
noSurfInds(data.elevation>-85)=0;
noSurfCols=find(noSurfInds==1);

for ii=1:length(noSurfCols)
    reflCol=reflTemp(:,noSurfCols(ii));
    nanCol=find(isnan(reflCol));
    diffCol=diff(nanCol);
    bigInd=min(find(diffCol>10));
    firstInd=nanCol(bigInd);
    maskTemp(1:firstInd,noSurfCols(ii))=3;
end
maskTemp=flipud(maskTemp);
echoMask(isnan(echoMask))=maskTemp(isnan(echoMask));
end
