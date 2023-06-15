% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function data=f_meltLayer_advanced(data,thresholds,figdir)

%debugFig=0;

%% Find zero degree altitude altitudes and indices

% disp('Searching zero degree altitude ...');
% 
% oneGate=data.range(2)-data.range(1);
% 
% [layerAlts,layerInds]=zeroDegIso(data);

%% Truncate to non missing and regions with sub 7 deg temps
gapSecs=10;
tempMax=7;
nonMissingInds=findNonMissingInds(data,gapSecs);

% Temperature too high
minTemp=min(data.TEMP,[],1,'omitnan');
nonMissingInds(minTemp>tempMax)=0;

dataInVars=fields(data);

dataShort=[];
for ii=1:length(dataInVars)
    dataShort.(dataInVars{ii})=data.(dataInVars{ii})(:,nonMissingInds==1);
end

%% Velocity derivative

% Remove upward motion to limit false detections
upInds=find(dataShort.elevation>=0);
downInds=find(dataShort.elevation<0);

velNoNeg=dataShort.VEL_MASKED;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);
velNoNeg(velNoNeg<=0)=nan;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);

% Median filter to smoth data
velNoNeg=medfilt2(velNoNeg,[3,7]);

% Derivative
velDiffShort=diff(velNoNeg,1);

% Add row of nans to make number of gates match with other variables
velDiffUp=velDiffShort(:,upInds);
velDiffUp=cat(1,nan(1,size(velDiffUp,2)),velDiffUp);
velDiffDown=velDiffShort(:,downInds);
velDiffDown=cat(1,velDiffDown,nan(1,size(velDiffDown,2)));
velDiff=nan(size(dataShort.VEL_MASKED));
velDiff(:,upInds)=velDiffUp;
velDiff(:,downInds)=velDiffDown;

%% Reflectivity derivative

% Median filter to smoth data
dataShort.DBZ_MASKED=medfilt2(dataShort.DBZ_MASKED,[3,7]);

% Erode in the vertical to get rid of edges
dbzMask=~isnan(dataShort.DBZ_MASKED);
dbzMask=imerode(dbzMask,strel('line',15,90));
dbzEroded=dataShort.DBZ_MASKED;
dbzEroded(dbzMask==0)=nan;

% Derivative
dbzDiffShort=diff(dbzEroded,1);

% Add row of nans to make number of gates match with other variables
dbzDiffUp=dbzDiffShort(:,upInds);
dbzDiffUp=cat(1,nan(1,size(dbzDiffUp,2)),dbzDiffUp);
dbzDiffDown=dbzDiffShort(:,downInds);
dbzDiffDown=cat(1,dbzDiffDown,nan(1,size(dbzDiffDown,2)));
dbzDiff=nan(size(dataShort.DBZ_MASKED));
dbzDiff(:,upInds)=dbzDiffUp;
dbzDiff(:,downInds)=dbzDiffDown;

% Reverse sign in upward
dbzDiff(:,upInds)=-dbzDiff(:,upInds);

%% Fuzzy logic to determine where melting layer is likely

% Get melting layer probability
meltProb=findMeltProb(dataShort,velDiff,dbzDiff);
meltProb(isnan(dataShort.DBZ_MASKED))=nan;
% Hard censor on temperature
meltProb(dataShort.TEMP<-1 | dataShort.TEMP>7)=nan;

%% Find areas where melting layer is very likely
maskHigh=meltProb>thresholds.meltProbHigh;
maskHigh=bwareaopen(maskHigh,35); % Remove small

%% Find high quality maxima
% Input fields for finding range gate with maximum
probHigh=meltProb;
probHigh(maskHigh==0)=nan;

% Maxima
[probMax,probMaxInd]=max(probHigh,[],1,'omitnan');

% Get altitudes of maxima
probMaxGetInd=probMaxInd(~isnan(probMax));
colprob=1:length(dataShort.time);
colprob(isnan(probMax))=[];
linprob=sub2ind(size(dataShort.VEL_MASKED),probMaxGetInd,colprob);
maxProbAlt=dataShort.asl(linprob);

%% Medians and stds of maxima to clean them up
smoothVal=201;
stdThresh=50;
maxProbAltNan=nan(size(dataShort.time));
maxProbAltNan(colprob)=maxProbAlt;

% Std
stdProb=movstd(maxProbAltNan,smoothVal,'omitnan');
stdProb(isnan(maxProbAltNan))=nan;
% Remove regions with high stds
maxProbAltNan(stdProb>stdThresh)=nan;

% Median
medProb=movmedian(maxProbAltNan,smoothVal,'omitnan');
medProb(isnan(maxProbAltNan))=nan;

%% Find high quality melt layer altitude
% Connect medians
medFilled=fillmissing(medProb,'linear');

% Keep only in close proximity
medMask=~isnan(medProb);
medMask=imdilate(medMask,strel('line',5000,0));
medFilled(medMask==0)=nan;

%% Melting layer mask, first guess
medFilledMat=repmat(medFilled,size(dataShort.DBZ_MASKED,1),1);
goodAlts=dataShort.asl-medFilledMat;
goodAlts(abs(goodAlts)>150)=nan;

meltMaskInit=(meltProb>thresholds.meltProbLow & ~isnan(goodAlts));
meltMaskInit=bwareaopen(meltMaskInit,35);

%% Find top and bottom altitudes
meltMaskTurned=meltMaskInit';
mFirst=size(meltMaskTurned,2);
[valFirst,locFirst]=max(meltMaskTurned,[],2);
kFirst=mFirst+1-locFirst;
kFirst(valFirst==0)=nan;
kFirst=abs(size(meltMaskInit,1)-kFirst);
linFirst=sub2ind(size(meltMaskInit),kFirst',1:size(meltMaskInit,2));
linFirst(isnan(linFirst))=[];
altFirstNoNan=dataShort.asl(linFirst);
altFirst=nan(size(dataShort.time));
altFirst(~isnan(kFirst))=altFirstNoNan;

mLast=size(meltMaskTurned,2);
[valLast,locLast]=max(fliplr(meltMaskTurned),[],2);
kLast=mLast+1-locLast;
kLast(valLast==0)=nan;
linLast=sub2ind(size(meltMaskInit),kLast',1:size(meltMaskInit,2));
linLast(isnan(linLast))=[];
altLastNoNan=dataShort.asl(linLast);
altLast=nan(size(dataShort.time));
altLast(~isnan(kLast))=altLastNoNan;

% Smooth slightly
smoothVal2=15;
altFirstS=movmedian(altFirst,smoothVal2,'omitnan');
altLastS=movmedian(altLast,smoothVal2,'omitnan');

% Remove edges
flMask=~isnan(altFirst);
flMask=imclose(flMask,strel('line',5,0));
altFirstS(flMask==0)=nan;
altLastS(flMask==0)=nan;

%% Final mask
maskInds=find(flMask==1);
meltMask=zeros(size(meltMaskInit));
for jj=1:length(maskInds)
    altCol=dataShort.asl(:,maskInds(jj));
    if dataShort.elevation(maskInds(jj))>=0
        rightInds=find(altCol<=altLastS(maskInds(jj)) & altCol>=altFirstS(maskInds(jj)));
    else
        rightInds=find(altCol>=altLastS(maskInds(jj)) & altCol<=altFirstS(maskInds(jj)));
    end
    meltMask(rightInds,maskInds(jj))=1;
end

%meltMask=imfill(meltMask,'holes');

%% Add censored data back in
data.velDiff=nan(size(data.DBZ_MASKED));
data.velDiff(:,nonMissingInds==1)=velDiff;

data.dbzDiff=nan(size(data.DBZ_MASKED));
data.dbzDiff(:,nonMissingInds==1)=dbzDiff;

data.meltMask=zeros(size(data.DBZ_MASKED));
data.meltMask(:,nonMissingInds==1)=meltMask;

data.meltProb=nan(size(data.DBZ_MASKED));
data.meltProb(:,nonMissingInds==1)=meltProb;

plotMedStd=0;
if plotMedStd
    ylimits=[0,8];
    meltTestPlot3;
end
end