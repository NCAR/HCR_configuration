% Flag melting layer
% 0=zero degree altitude
% 1=melting layer detected
% 2=melting layer interpolated
% 3=melting layer defined as zero degree altitude
function [BBfinishedOut iceLevAsl offset]= f_meltLayer_advanced(data,zeroAdjustMeters,thresholds,figdir)

debugFig=0;

%% Find zero degree altitude altitudes and indices

disp('Searching zero degree altitude ...');

oneGate=data.range(2)-data.range(1);

[layerAlts,layerInds]=zeroDegIso(data);

%% Truncate to non missing and regions with sub deg temps
gapSecs=10;
tempMinMax=[-1,7];
nonMissingInds=findNonMissingInds(data,gapSecs);

% Temperature too high
minTemp=min(data.TEMP,[],1,'omitnan');
nonMissingInds(minTemp>tempMinMax(2))=0;

dataInVars=fields(data);

dataShort=[];
for ii=1:length(dataInVars)
    dataShort.(dataInVars{ii})=data.(dataInVars{ii})(:,nonMissingInds==1);
end

%% Prepare VEL data

% Median filter
dataShort.VEL_MASKED=medfilt2(dataShort.VEL_MASKED,[3,7]);

upInds=find(dataShort.elevation>=0);
downInds=find(dataShort.elevation<0);

velNoNeg=dataShort.VEL_MASKED;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);
velNoNeg(velNoNeg<=0)=nan;
velNoNeg(:,upInds)=-velNoNeg(:,upInds);

%% Velocity derivative
velDiffShort=diff(velNoNeg,1);
velDiffUp=velDiffShort(:,upInds);
velDiffUp=cat(1,nan(1,size(velDiffUp,2)),velDiffUp);
velDiffDown=velDiffShort(:,downInds);
velDiffDown=cat(1,velDiffDown,nan(1,size(velDiffDown,2)));
velDiff=nan(size(dataShort.VEL_MASKED));
velDiff(:,upInds)=velDiffUp;
velDiff(:,downInds)=velDiffDown;

%% Fuzzy logic to determine where melting layer is likely

meltProb=findMeltProb(dataShort,velDiff);
meltProbThreshHigh=0.7;
meltProbThreshLow=0.55;
meltProb(isnan(dataShort.DBZ_MASKED))=nan;

maskHigh=meltProb>meltProbThreshHigh;
maskHigh=bwareaopen(maskHigh,35);

%% Input fields
velField=velDiff;
velField(maskHigh==0)=nan;
ldrField=dataShort.LDR;
ldrField(maskHigh==0)=nan;

%% Maxima
[velMax,velMaxInd]=max(velField,[],1,'omitnan');
[ldrMax,ldrMaxInd]=max(ldrField,[],1,'omitnan');

%% Get altitudes of max
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

%% Medians and stds
smoothVal=201;
maxVelAltNan=nan(size(dataShort.time));
maxVelAltNan(colvel)=maxVelAlt;
% Std
stdVel=movstd(maxVelAltNan,smoothVal,'omitnan');
stdVel(isnan(maxVelAltNan))=nan;
maxVelAltNan(stdVel>100)=nan;
% Median
medVel=movmedian(maxVelAltNan,smoothVal,'omitnan');
medVel(isnan(maxVelAltNan))=nan;

maxLdrAltNan=nan(size(dataShort.time));
maxLdrAltNan(colldr)=maxLdrAlt;
% Std
stdLdr=movstd(maxLdrAltNan,smoothVal,'omitnan');
stdLdr(isnan(maxLdrAltNan))=nan;
maxLdrAltNan(stdLdr>100)=nan;
% Median
medLdr=movmedian(maxLdrAltNan,smoothVal,'omitnan');
medLdr(isnan(maxLdrAltNan))=nan;
% Combine median and connect
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

meltMask=(meltProb>meltProbThreshLow & ~isnan(goodAlts));
meltMask=imclose(meltMask,strel('disk',15));
meltMask=imfill(meltMask,'holes');
meltMask=bwareaopen(meltMask,35);

%% Plot
disp('Plotting ...')
ylimits=[0,5];

showPlot='off';

close all

meltTestPlot1;
meltTestPlot2;

end