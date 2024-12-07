function [echoMask,antStat] = echoMask(data,ls)

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

echoMask=nan(size(data.(['DBZ',ls])));

%% Find bang (6)
if strcmp(ls,'_short')
    echoMask(1:17,:)=6;
else
    echoMask(1:19,:)=6;
end

%% Missing data (11)
bangData=data.(['DBZ',ls])(1:14,:);
bangData(bangData<-10)=nan;
bangMask=zeros(size(bangData));
bangMask(isnan(bangData))=1;
sumNan=sum(bangMask,1,'omitnan');

missInds=find(sumNan>13);
echoMask(:,missInds)=11;

%% NS cal (10)
firstGate=data.(['DBMVC',ls])(1,:);
nscalInds=find(firstGate>-96 & firstGate<-87);
nscalMask=zeros(size(firstGate));
nscalMask(nscalInds)=1;

echoMask(:,find(nscalMask==1))=10;

%% Cloud echo (1)

refl=data.(['DBZ',ls]);
refl(echoMask>0)=nan;

goodInds=find(~isnan(refl));
echoMask(goodInds)=1;

%% Speckle (2)
specCut=100;

reflTemp=data.(['DBZ',ls]);
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
end
