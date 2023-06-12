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

%% Truncate to non missing and regions with sub deg temps
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

%% Prepare VEL data

% Median filter to smoth data
dataShort.VEL_MASKED=medfilt2(dataShort.VEL_MASKED,[3,7]);

% Remove upward motion to limit false detections
upInds=find(dataShort.elevation>=0);
downInds=find(dataShort.elevation<0);

velNoNeg=dataShort.VEL_MASKED;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);
velNoNeg(velNoNeg<=0)=nan;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);

%% Velocity derivative
velDiffShort=diff(velNoNeg,1);

% Add row of nans to make number of gates match with other variables
velDiffUp=velDiffShort(:,upInds);
velDiffUp=cat(1,nan(1,size(velDiffUp,2)),velDiffUp);
velDiffDown=velDiffShort(:,downInds);
velDiffDown=cat(1,velDiffDown,nan(1,size(velDiffDown,2)));
velDiff=nan(size(dataShort.VEL_MASKED));
velDiff(:,upInds)=velDiffUp;
velDiff(:,downInds)=velDiffDown;

%% Fuzzy logic to determine where melting layer is likely

% Get melting layer probability
meltProb=findMeltProb(dataShort,velDiff);
meltProb(isnan(dataShort.DBZ_MASKED))=nan;

% Threshold for almost certain
thresholds.meltProbHigh=0.7;
% Threshold to fill in areas
thresholds.meltProbLow=0.55;

%% Find areas where melting layer is very likely
maskHigh=meltProb>thresholds.meltProbHigh;
maskHigh=bwareaopen(maskHigh,35); % Remove small

%% Find high quality maxima
% Input fields for finding range gate with maximum
velField=velDiff;
velField(maskHigh==0)=nan;
ldrField=dataShort.LDR;
ldrField(maskHigh==0)=nan;

% Maxima
[velMax,velMaxInd]=max(velField,[],1,'omitnan');
[ldrMax,ldrMaxInd]=max(ldrField,[],1,'omitnan');

% Get altitudes of maxima
velMaxPlotInd=velMaxInd(~isnan(velMax));
colvel=1:length(dataShort.time);
colvel(isnan(velMax))=[];
linvel=sub2ind(size(dataShort.VEL_MASKED),velMaxPlotInd,colvel);
maxVelAlt=dataShort.asl(linvel);

ldrMaxPlotInd=ldrMaxInd(~isnan(ldrMax));
colldr=1:length(dataShort.time);
colldr(isnan(ldrMax))=[];
linldr=sub2ind(size(dataShort.LDR),ldrMaxPlotInd,colldr);
maxLdrAlt=dataShort.asl(linldr);

%% Medians and stds of maxima to clean them up
smoothVal=201;
stdThresh=100;
maxVelAltNan=nan(size(dataShort.time));
maxVelAltNan(colvel)=maxVelAlt;

% Std vel
stdVel=movstd(maxVelAltNan,smoothVal,'omitnan');
stdVel(isnan(maxVelAltNan))=nan;
% Remove regions with high stds
maxVelAltNan(stdVel>stdThresh)=nan;

% Median vel
medVel=movmedian(maxVelAltNan,smoothVal,'omitnan');
medVel(isnan(maxVelAltNan))=nan;

maxLdrAltNan=nan(size(dataShort.time));
maxLdrAltNan(colldr)=maxLdrAlt;

% Std ldr
stdLdr=movstd(maxLdrAltNan,smoothVal,'omitnan');
stdLdr(isnan(maxLdrAltNan))=nan;
% Remove regions with high stds
maxLdrAltNan(stdLdr>stdThresh)=nan;

% Median ldr
medLdr=movmedian(maxLdrAltNan,smoothVal,'omitnan');
medLdr(isnan(maxLdrAltNan))=nan;

%% Find high quality melt layer altitude
% Combine vel and ldr medians and connect
medComb=mean([medVel;medLdr],1,'omitnan');
medFilled=fillmissing(medComb,'linear');

% Keep only in close proximity
medMask=~isnan(medComb);
medMask=imdilate(medMask,strel('line',5000,0));
medFilled(medMask==0)=nan;

%% Melting layer mask
medFilledMat=repmat(medFilled,size(dataShort.DBZ_MASKED,1),1);
goodAlts=dataShort.asl-medFilledMat;
goodAlts(abs(goodAlts)>150)=nan;

meltMask=(meltProb>thresholds.meltProbLow & ~isnan(goodAlts));
meltMask=bwareaopen(meltMask,35);
meltMask=imfill(meltMask,'holes');

%% Add censored data back in
data.velDiff=nan(size(data.DBZ_MASKED));
data.velDiff(:,nonMissingInds==1)=velDiff;

data.meltMask=zeros(size(data.DBZ_MASKED));
data.meltMask(:,nonMissingInds==1)=meltMask;

data.meltProb=nan(size(data.DBZ_MASKED));
data.meltProb(:,nonMissingInds==1)=meltProb;

plotMedStd=1;
if plotMedStd
    ylimits=[0,8];
    meltTestPlot3;
end
end