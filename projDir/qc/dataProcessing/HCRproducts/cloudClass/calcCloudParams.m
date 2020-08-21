function cloudParams=calcCloudParams(dataCut);
% Calculate parameters for all individual clouds

cloudParams=[];

%% Max and min altitude, temperature, and precipitation

% Convert from ASL to AGL
agl=dataCut.asl-dataCut.TOPO./1000;

% Mean temperature
%cloudParams.meanTemp=mean(dataCut.TEMP,'omitnan');

% Initialize min/max agl, temp, and precip
minAglAll=nan(1,size(dataCut.cloudPuzzle,2));
maxAglAll=nan(1,size(dataCut.cloudPuzzle,2));
minTempAll=nan(1,size(dataCut.cloudPuzzle,2));
maxTempAll=nan(1,size(dataCut.cloudPuzzle,2));
precip=nan(1,size(dataCut.cloudPuzzle,2));

for jj=1:size(dataCut.cloudPuzzle,2)
    % Min/max alt
    aglCol=agl(:,jj);
    [minAglAll(jj) minIndagl]=min(aglCol);
    [maxAglAll(jj) maxIndagl]=max(aglCol);
    
    % Temperature
    minTempAll(jj)=dataCut.TEMP(minIndagl,jj);
    maxTempAll(jj)=dataCut.TEMP(maxIndagl,jj);
    
    % Check if flying in cloud
    % Pointing down
    if dataCut.elevation(jj)<0 & maxAglAll(jj)<10000 & maxIndagl==18
        maxAglAll(jj)=nan;
        maxTempAll(jj)=nan;
    end
    % Pointing up
    if dataCut.elevation(jj)>=0 & minAglAll(jj)<10000 & minIndagl==18
        minAglAll(jj)=nan;
        minTempAll(jj)=nan;
    end
    
    % Precip
    if dataCut.elevation(jj)<0
        surfInd=min(find(dataCut.FLAG(:,jj)==7));
        if ~isempty(surfInd)
            precip(jj)=mean(dataCut.puzzleOne(surfInd-5:surfInd-1,jj));
        elseif any(dataCut.flagCensored(:,jj)==3)
            precip(jj)=1;
        end
    end
end

percWanted=0.02;

% AGL: Make sure we have enough data and calculate percentiles
if length(find(~isnan(minAglAll)))>length(minAglAll)/2
    sortedMin=sort(minAglAll,'ascend');
    percIndMin=round(percWanted*length(minAglAll));
    cloudParams.minAgl=sortedMin(percIndMin);
    cloudParams.meanMinAgl=mean(minAglAll,'omitnan');
    cloudParams.stdMinAgl=std(minAglAll,'omitnan');
else
    cloudParams.minAgl=nan;
    cloudParams.meanMinAgl=nan;
    cloudParams.stdMinAgl=nan;
end

if length(find(~isnan(maxAglAll)))>length(maxAglAll)/2
    sortedMax=sort(maxAglAll,'descend');
    percIndMax=round(percWanted*length(maxAglAll));
    cloudParams.maxAgl=sortedMax(percIndMax);
    cloudParams.meanMaxAgl=mean(maxAglAll,'omitnan');
    cloudParams.stdMaxAgl=std(maxAglAll,'omitnan');
else
    cloudParams.maxAgl=nan;
    cloudParams.meanMaxAgl=nan;
    cloudParams.stdMaxAgl=nan;
end

% TEMP: Make sure we have enough data and calculate percentiles
if length(find(~isnan(minTempAll)))>length(minTempAll)/2
    sortedMinTemp=sort(minTempAll,'ascend');
    percIndMinTemp=round(percWanted*length(minTempAll));
    %cloudParams.minBaseTemp=sortedMinTemp(percIndMinTemp);
    cloudParams.meanBaseTemp=mean(minTempAll,'omitnan');
else
    %cloudParams.minBaseTemp=nan;
    cloudParams.meanBaseTemp=nan;
end

if length(find(~isnan(maxTempAll)))>length(maxTempAll)/2
    sortedMaxTemp=sort(maxTempAll,'ascend');
    percIndMaxTemp=round(percWanted*length(maxTempAll));
    cloudParams.minTopTemp=sortedMaxTemp(percIndMaxTemp);
    cloudParams.meanTopTemp=mean(maxTempAll,'omitnan');
else
    cloudParams.minTopTemp=nan;
    cloudParams.meanTopTemp=nan;
end

% Precip
cloudParams.numPrecip=length(find(~isnan(precip)));
if cloudParams.numPrecip>length(precip)*0.03
    precCloud=1;
else
    precCloud=0;
end

cloudParams.precip=precCloud;

% Mean max refl and mean max refl height, mean temp at max refl height
[maxRefl maxReflInds]=max(dataCut.DBZ,[],1);
cloudParams.meanMaxRefl=mean(maxRefl,'omitnan');
cloudParams.stdMaxRefl=std(maxRefl,'omitnan');

% Total maximum refl
sortedMaxRefl=sort(maxRefl,'descend');
percIndMaxRefl=round(percWanted*length(maxRefl));
cloudParams.maxMaxRefl=sortedMaxRefl(percIndMaxRefl);

% Temperature and height of max refl
wholeMaxReflInds=sub2ind(size(dataCut.DBZ),maxReflInds,1:size(dataCut.DBZ,2));
maxReflAgl=agl(wholeMaxReflInds);
maxReflTemp=dataCut.TEMP(wholeMaxReflInds);

cloudParams.meanMaxReflAgl=mean(maxReflAgl,'omitnan');
cloudParams.meanMaxReflTemp=mean(maxReflTemp,'omitnan');

% Mean latitude
cloudParams.meanLat=mean(dataCut.latitude,'omitnan');

% Cloud length in km (make sure length is not cut off because of a/descent
% or missing data
sumDBZ=sum(dataCut.DBZ,1,'omitnan');
startColZ=dataCut.DBZ(:,2);
startColF=dataCut.FLAG(:,1);
startColF(isnan(startColZ))=nan;

endColZ=dataCut.DBZ(:,end-1);
endColF=dataCut.FLAG(:,end);
endColF(isnan(endColZ))=nan;

if sumDBZ(1)~=0 | sumDBZ(end)~=0 | sum(startColF,'omitnan')~=0 | sum(endColF,'omitnan')~=0
    cloudParams.lengthKM=nan;
else
    [cloudParams.lengthKM ~]=lldistkm([dataCut.latitude(1) dataCut.longitude(1)],[dataCut.latitude(end) dataCut.longitude(end)]);
end

% Maximum 10 dBZ height
tenDBZ=dataCut.DBZ;
tenDBZ(tenDBZ<10)=nan;

agl10dbz=agl(~isnan(tenDBZ));
if isempty(agl10dbz)
    cloudParams.max10dbzAgl=nan;
else
    sortedAgl10=sort(agl10dbz,'descend');
    percIndAgl10=round(percWanted*length(sortedAgl10));
    cloudParams.max10dbzAgl=sortedAgl10(percIndAgl10);
end

% Cloud thickness
cloudThick=abs(maxAglAll-minAglAll);
cloudParams.meanThickness=mean(cloudThick,'omitnan');

sortedThick=sort(cloudThick,'descend');
percIndThick=round(percWanted*length(sortedThick));
cloudParams.maxThickness=sortedThick(percIndThick);

% Inhomogeneity
maxReflLin=10.^(maxRefl./10);
cloudParams.inhomo=std(maxReflLin,'omitnan')/mean(maxReflLin,'omitnan');

end