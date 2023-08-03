% Flag melting layer
% 0=cold
% 1=melting
% 2=warm
function data=f_meltLayer_advanced(data,thresholds,figdir)

disp('Preparing fields ...')

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
upIndsS=find(dataShort.elevation>=0);
downIndsS=find(dataShort.elevation<0);

velNoNeg=dataShort.VEL_MASKED;
velNoNeg(:,upIndsS)=-velNoNeg(:,upIndsS);
velNoNeg(velNoNeg<=0)=nan;
velNoNeg(:,upIndsS)=-velNoNeg(:,upIndsS);

% Median filter to smoth data
velNoNeg=medfilt2(velNoNeg,[3,7]);

% Derivative
velDiffShort=diff(velNoNeg,1);

% Add row of nans to make number of gates match with other variables
velDiffUp=velDiffShort(:,upIndsS);
velDiffUp=cat(1,nan(1,size(velDiffUp,2)),velDiffUp);
velDiffDown=velDiffShort(:,downIndsS);
velDiffDown=cat(1,velDiffDown,nan(1,size(velDiffDown,2)));
velDiff=nan(size(dataShort.VEL_MASKED));
velDiff(:,upIndsS)=velDiffUp;
velDiff(:,downIndsS)=velDiffDown;

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
dbzDiffUp=dbzDiffShort(:,upIndsS);
dbzDiffUp=cat(1,nan(1,size(dbzDiffUp,2)),dbzDiffUp);
dbzDiffDown=dbzDiffShort(:,downIndsS);
dbzDiffDown=cat(1,dbzDiffDown,nan(1,size(dbzDiffDown,2)));
dbzDiff=nan(size(dataShort.DBZ_MASKED));
dbzDiff(:,upIndsS)=dbzDiffUp;
dbzDiff(:,downIndsS)=dbzDiffDown;

% Reverse sign in upward
dbzDiff(:,upIndsS)=-dbzDiff(:,upIndsS);

%% Fuzzy logic to determine where melting layer is likely

disp('Fuzzy logic ...')

% Get melting layer probability
meltProb=findMeltProb(dataShort,velDiff,dbzDiff);
meltProb(isnan(dataShort.DBZ_MASKED))=nan;
% Hard censor on temperature
meltProb(dataShort.TEMP<-1 | dataShort.TEMP>7)=nan;

% Remove small temp strips
tempMaskIn=dataShort.TEMP>0;
tempMask=imclose(tempMaskIn,strel('line',7,90));
tempMask=~(bwareaopen(~tempMask,50000));
tempMask=double(tempMask);
meltProb(tempMaskIn==0 & tempMask==1)=nan;

%% Find areas where melting layer is very likely
maskHigh=meltProb>thresholds.meltProbHigh;
maskHigh=bwareaopen(maskHigh,40); % Remove small

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
maxProbAltNan(abs(medProb-maxProbAltNan)>1000)=nan;
medProb(isnan(maxProbAltNan))=nan;

%% Find high quality melt layer altitude
% Connect medians
medFilled=fillmissing(medProb,'linear');

% Keep only in close proximity
medMask=~isnan(medProb);
medMask=imdilate(medMask,strel('line',5000,0));
medFilled(medMask==0)=nan;

%% Melting layer mask
medFilledMat=repmat(medFilled,size(dataShort.DBZ_MASKED,1),1);
goodAlts=dataShort.asl-medFilledMat;
goodAlts(abs(goodAlts)>150)=nan;

meltMask=(meltProb>thresholds.meltProbLow & ~isnan(goodAlts));

% Remove small over extinct
extinctMask=any(dataShort.FLAG==3,1);
extinctMask=imclose(extinctMask,strel('line',15,0));
extinctMask=imdilate(extinctMask,strel('line',55,0));

meltMask1D=any(meltMask==1,1);
meltMask1D=imclose(meltMask1D,strel('line',55,0));
meltMask1D=bwareaopen(meltMask1D,500);

meltMask(:,(extinctMask==1 & meltMask1D==0))=0;

% Clean up
meltMask=bwareaopen(meltMask,35); % Remove small
meltMask=imclose(meltMask,strel('disk',25)); % Smooth and join
meltMask=imfill(meltMask,'holes'); % Remove holes

%% Interpolate melting layer altitude
% Get max probability within mask
probMask=meltProb;
probMask(meltMask==0)=nan;
noNan=any(~isnan(probMask),1);
[~,maxInd2]=max(probMask,[],1,'omitnan');
maxInd2lin=sub2ind(size(probMask),maxInd2,1:size(probMask,2));
maxInd2lin(noNan==0)=[];

% Get altitude of maximum
maxAlt2noNan=dataShort.asl(maxInd2lin);
maxAlt2=nan(size(dataShort.time));
maxAlt2(noNan==1)=maxAlt2noNan;

% Smooth slightly
maxAltS=movmedian(maxAlt2,15,'omitnan');
maxAltS(isnan(maxAlt2))=nan;

% Connect
maxAltIntShort=fillmissing(maxAltS,'linear','EndValues','nearest');
% Only in close proximity
transLength=1800; % 3 minutes (1800 1/10 seconds)
largerNoNan=imdilate(noNan,strel('line',transLength,0)); % 1.5 minutes, imdilate makes it half
% Interpolate over holes that are twice the transition length
largerNoNan=imclose(largerNoNan,strel('line',round(transLength+transLength/2),0)); % 3 minutes, imclose takes whole, plus one 1.5 mins
maxAltIntShort(largerNoNan==0)=nan;

%% Add censored data back in
data.velDiff=nan(size(data.DBZ_MASKED));
data.velDiff(:,nonMissingInds==1)=velDiff;

data.dbzDiff=nan(size(data.DBZ_MASKED));
data.dbzDiff(:,nonMissingInds==1)=dbzDiff;

data.meltMask=zeros(size(data.DBZ_MASKED));
data.meltMask(:,nonMissingInds==1)=meltMask;

data.meltProb=nan(size(data.DBZ_MASKED));
data.meltProb(:,nonMissingInds==1)=meltProb;

maxAltInt=nan(size(data.time));
maxAltInt(nonMissingInds==1)=maxAltIntShort;

maxAltSall=nan(size(data.time));
maxAltSall(nonMissingInds==1)=maxAltS;

%% Find warm and cold regions

disp('Finding warm and cold regions of troposphere ...')
[warmCold,zeroAltsAdj,meltOnly]=findWarmCold(data);

%% Merge interpolated melt alt and adjusted freezing levels

disp('Matching found melting layers with zero deg layers ...')

% Only keep zero alts that go from cold to warm
zeroAltMelt=zeroAltsAdj;
zeroAltMelt(meltOnly==0,:)=nan;

% Handle max alt prep
% Find connected
maxAltIntMask=~isnan(maxAltInt);
maxAltRegAll=bwconncomp(maxAltIntMask);

% Initiate output
zeroAltReal=zeroAltsAdj; % Zero deg

% Loop through individual found melting layers, put zero deg and found
% melting layer together with some distance in-between
rightIndAll=[];
zeroAltAdd=nan(size(zeroAltReal));
for jj=1:maxAltRegAll.NumObjects
    maxAltReg=maxAltInt(maxAltRegAll.PixelIdxList{jj}); % Found altitude
    zeroMeltReg=zeroAltMelt(:,maxAltRegAll.PixelIdxList{jj}); % Zero deg alt
    % Horizontal median to find the correct zero deg alt index
    maxMed=median(maxAltReg,'omitnan');
    zeroMed=median(zeroMeltReg,2,'omitnan');
    diffAlt=zeroMed-maxMed;
    % Find index
    [~,rightInd]=min(abs(diffAlt));
    if ~isnan(rightInd)
        rightIndAll=[rightIndAll;rightInd];
    end
    % Make room for detected
    minPixReg=max([1,min(maxAltRegAll.PixelIdxList{jj})-transLength]);
    maxPixReg=min([max(maxAltRegAll.PixelIdxList{jj})+transLength,size(zeroAltReal,2)]);
    zeroAltReal(rightInd,minPixReg:maxPixReg)=nan;
    % Add found melting layer
    zeroAltAdd(rightInd,maxAltRegAll.PixelIdxList{jj})=maxAltReg;
end

zeroAltReal(~isnan(zeroAltAdd))=zeroAltAdd(~isnan(zeroAltAdd));

% Combine
for kk=1:length(rightIndAll)
    % Close wholes
    thisAlt=zeroAltReal(rightIndAll(kk),:);
    thisMask=~isnan(thisAlt);
    thisMask=imdilate(thisMask,strel('line',transLength,0));
    thisAlt=fillmissing(thisAlt,'linear','EndValues','nearest');
    % But only in vicinity
    thisAlt(thisMask==0)=nan;
    zeroAltReal(rightIndAll(kk),:)=thisAlt;
end

% Extend in climbs and descents
for kk=1:size(zeroAltReal,1)
    thisAlt=zeroAltReal(kk,:);
    thisMask=~isnan(thisAlt);
    thisMask=imdilate(thisMask,strel('line',500,0));
    thisAltFilled=fillmissing(thisAlt,'linear','EndValues','nearest');
    thisAltFilled(thisMask==0)=nan;
    zeroAltReal(kk,:)=thisAltFilled;
end

plotYes=0;
if plotYes
    for ii=1:size(zeroAltsAdj,1)
        plot(data.time,zeroAltReal(ii,:)./1000,'-m','LineWidth',2);
    end
    ax=gca;
    ax.SortMethod='childorder';
end

%% Create warm/cold/melting mat

disp('Creating melting layer output ...')

meltLayerOut=warmCold;
for ll=1:length(data.time)
    colAlts=[zeroAltReal(:,ll),meltOnly];
    colAlts(isnan(colAlts(:,1)),:)=[];
    if isempty(colAlts)
        continue
    end
    altCol=data.asl(:,ll);
    meltCol=nan(size(altCol));
    colAltsSort=sortrows(colAlts);
    colAltsSort=cat(1,[0,0],colAltsSort);
    colAltsSort=cat(1,colAltsSort,[20000,abs(colAltsSort(end,2)-1)]);
    for mm=1:size(colAltsSort,1)-1
        if colAltsSort(mm+1,2)==1
            meltCol(altCol>=colAltsSort(mm,1) & altCol<=colAltsSort(mm+1,1))=2;
        elseif colAltsSort(mm+1,2)==0
            meltCol(altCol>=colAltsSort(mm,1) & altCol<=colAltsSort(mm+1,1))=0;
        end
        meltLayerOut(:,ll)=meltCol;
    end
end

% Sanity check
meltLayerOut(data.TEMP<-1)=0;
meltLayerOut(data.TEMP>7)=2;

% Clean up small areas
coldMask=meltLayerOut==0;
coldMask=bwareaopen(coldMask,200000);
coldMask=imclose(coldMask,strel('disk',1));
meltLayerOut(coldMask==0)=2;
meltLayerOut(coldMask==1)=0;

warmMask=meltLayerOut==2;
warmMask=bwareaopen(warmMask,200000);
warmMask=imclose(warmMask,strel('disk',1));
meltLayerOut(warmMask==0)=0;
meltLayerOut(warmMask==1)=2;

% Final numbers
meltLayerOut(meltLayerOut==0)=21; % Cold
meltLayerOut(meltLayerOut==2)=9; % Warm

% Split melt mask into warm and cold
meltMaskAlt=data.asl;
meltMaskAlt(data.meltMask==0)=nan;

maxMat=repmat(maxAltSall,size(data.DBZ_MASKED,1),1);
meltMaskDiff=meltMaskAlt-maxMat;

meltMaskSplit=nan(size(data.meltMask));
meltMaskSplit(meltMaskDiff>=0)=19; % Melting cold
meltMaskSplit(meltMaskDiff<0)=11; % Melting warm

% Add melting areas
meltLayerOut(~isnan(meltMaskSplit))=meltMaskSplit(~isnan(meltMaskSplit));
data.meltLayer=meltLayerOut;

data.meltLayer(isnan(data.TEMP))=nan;
data.meltLayer(data.FLAG~=1)=nan;

%% Testing plot for stds and medians
plotMedStd=0;
if plotMedStd
    ylimits=[0,8];
    meltTestPlot3;
end
end