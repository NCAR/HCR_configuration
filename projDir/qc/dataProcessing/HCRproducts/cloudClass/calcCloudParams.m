function cloudParams=calcCloudParams(dataCut);
% Calculate parameters for all individual clouds

cloudParams=[];

debugFig=0;

if debugFig
    close all
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,900,800]);
    
    hold on;
    sub1=surf(dataCut.time,dataCut.asl,dataCut.DBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    %ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([dataCut.time(1),dataCut.time(end)]);
    title('Reflectivity')
    grid on
end

%% Max and min altitude, temperature, and precipitation

% Convert from ASL to AGL
agl=dataCut.asl-dataCut.TOPO./1000;

planeAlt=dataCut.altitude-dataCut.TOPO;

% Mean temperature
%cloudParams.meanTemp=mean(dataCut.TEMP,'omitnan');

% Initialize min/max agl, temp, and precip
minAglAll=nan(1,size(dataCut.cloudPuzzle,2));
maxAglAll=nan(1,size(dataCut.cloudPuzzle,2));
minTempAll=nan(1,size(dataCut.cloudPuzzle,2));
maxTempAll=nan(1,size(dataCut.cloudPuzzle,2));
minAglAllOrig=nan(1,size(dataCut.cloudPuzzle,2));
maxAglAllOrig=nan(1,size(dataCut.cloudPuzzle,2));
minTempAllOrig=nan(1,size(dataCut.cloudPuzzle,2));
maxTempAllOrig=nan(1,size(dataCut.cloudPuzzle,2));
precipIn=nan(1,size(dataCut.cloudPuzzle,2));

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
    maxAglAllOrig(jj)=maxAglAll(jj);
    maxTempAllOrig(jj)=maxTempAll(jj);
    if dataCut.elevation(jj)<0 & maxAglAll(jj)<10000 & maxIndagl==18
        maxAglAll(jj)=nan;
        maxTempAll(jj)=nan;
    end
    % Pointing up
    minAglAllOrig(jj)=minAglAll(jj);
    minTempAllOrig(jj)=minTempAll(jj);
    if dataCut.elevation(jj)>=0 & minAglAll(jj)<10000 & minIndagl==18 & abs(planeAlt(jj))>20
        minAglAll(jj)=nan;
        minTempAll(jj)=nan;
    end
    
    % Precip
    if ~isnan(minAglAll(jj)) | (abs(planeAlt(jj))<20 & dataCut.elevation>=0)
        if abs(planeAlt(jj))<20 & dataCut.elevation>=0
            surfInd=18;
        else
            surfInd=min(find(dataCut.FLAG(:,jj)==7));
        end
        if ~isempty(surfInd) & dataCut.elevation<0
            precipIn(jj)=mean(dataCut.puzzleOne(surfInd-5:surfInd-1,jj));
        elseif ~isempty(surfInd) & dataCut.elevation>=0
            precipIn(jj)=mean(dataCut.puzzleOne(surfInd+1:surfInd+5,jj));
        elseif any(dataCut.flagCensored(:,jj)==3)
            precipIn(jj)=1;
        end
    else
        precipIn(jj)=nan;
    end
end

percWanted=0.02;

% AGL: Make sure we have enough data and calculate percentiles
if length(find(~isnan(minAglAll)))>length(minAglAll)/2
    sortedMin=sort(minAglAllOrig,'ascend');
    sortedMin(isnan(sortedMin))=[];
    percIndMin=round(percWanted*length(minAglAllOrig));
    cloudParams.minAgl=sortedMin(percIndMin);
    cloudParams.meanMinAgl=mean(minAglAllOrig,'omitnan');
    cloudParams.stdMinAgl=std(minAglAll,'omitnan');
else
    cloudParams.minAgl=nan;
    cloudParams.meanMinAgl=nan;
    cloudParams.stdMinAgl=nan;
end

if length(find(~isnan(maxAglAll)))>length(maxAglAll)/2
    sortedMax=sort(maxAglAllOrig,'descend');
    sortedMax(isnan(sortedMax))=[];
    percIndMax=round(percWanted*length(maxAglAllOrig));
    cloudParams.maxAgl=sortedMax(percIndMax);
    cloudParams.meanMaxAgl=mean(maxAglAllOrig,'omitnan');
    cloudParams.stdMaxAgl=std(maxAglAll,'omitnan');
else
    cloudParams.maxAgl=nan;
    cloudParams.meanMaxAgl=nan;
    cloudParams.stdMaxAgl=nan;
end

% TEMP: Make sure we have enough data and calculate percentiles
if length(find(~isnan(minTempAll)))>length(minTempAll)/2
    sortedMinTemp=sort(minTempAllOrig,'ascend');
    sortedMinTemp(isnan(sortedMinTemp))=[];
    percIndMinTemp=round(percWanted*length(minTempAllOrig));
    %cloudParams.minBaseTemp=sortedMinTemp(percIndMinTemp);
    cloudParams.meanBaseTemp=mean(minTempAllOrig,'omitnan');
else
    %cloudParams.minBaseTemp=nan;
    cloudParams.meanBaseTemp=nan;
end

if length(find(~isnan(maxTempAll)))>length(maxTempAll)/2
    sortedMaxTemp=sort(maxTempAllOrig,'ascend');
    sortedMaxTemp(isnan(sortedMaxTemp))=[];
    percIndMaxTemp=round(percWanted*length(maxTempAllOrig));
    cloudParams.minTopTemp=sortedMaxTemp(percIndMaxTemp);
    cloudParams.meanTopTemp=mean(maxTempAllOrig,'omitnan');
else
    cloudParams.minTopTemp=nan;
    cloudParams.meanTopTemp=nan;
end

% Precip
cloudParams.numPrecip=length(find(~isnan(precipIn)));
if cloudParams.numPrecip>length(precipIn)*0.03
    precCloud=1;
elseif mean(dataCut.elevation<0) | cloudParams.minAgl>0.5
    precCloud=0;
else
    precCloud=nan;
end

cloudParams.precip=precCloud;

% Intense and very intense precip
cloudParams.intPrecip=0;
cloudParams.veryIntPrecip=0;

if cloudParams.precip==1
    oceanLand=nan(size(dataCut.time)); % Ocean=1, land=2, extinct=3
    oceanLand(any(dataCut.FLAG==7,1))=1;
    oceanLand(any(dataCut.FLAG==8,1))=2;
    
    oceanLand=fillmissing(oceanLand,'nearest');
        
    surfReflMasked=dataCut.surfRefl;
    surfReflMasked(isnan(precipIn))=nan;
    
    % Missing surface echo
    cloudParams.veryIntPrecip=cloudParams.veryIntPrecip+length(find(~isnan(precipIn) & isnan(surfReflMasked)));
    % Very intense
    cloudParams.veryIntPrecip=cloudParams.veryIntPrecip+length(find(surfReflMasked<-10));
    % Ocean
    cloudParams.intPrecip=cloudParams.intPrecip+length(find(oceanLand==1 & surfReflMasked<20));
    % Land
    cloudParams.intPrecip=cloudParams.intPrecip+length(find(oceanLand==2 & surfReflMasked<10));
end

% Mean max refl and mean max refl height, mean temp at max refl height
[maxRefl maxReflInds]=max(dataCut.DBZ,[],1);
maxRefl(isnan(maxRefl))=[];
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
testLength=dataCut.FLAG;
testLength(~isnan(dataCut.cloudPuzzle))=nan;
testLength(testLength==4)=nan;

sumDBZ=sum(dataCut.DBZ,1,'omitnan');
startColZ=dataCut.DBZ(:,2);
startColF=testLength(:,1);
startColF(isnan(startColZ))=nan;

endColZ=dataCut.DBZ(:,end-1);
endColF=testLength(:,end);
endColF(isnan(endColZ))=nan;

[cloudParams.lengthKM ~]=lldistkm([dataCut.latitude(1) dataCut.longitude(1)],[dataCut.latitude(end) dataCut.longitude(end)]);

if (sumDBZ(1)~=0 | sumDBZ(end)~=0 | sum(startColF,'omitnan')~=0 | sum(endColF,'omitnan')~=0) & cloudParams.lengthKM<=271
    cloudParams.lengthKM=nan;
end

% Maximum 10 dBZ height
tenDBZ=dataCut.DBZ;
tenDBZ(tenDBZ<10)=nan;

agl10dbz=agl(~isnan(tenDBZ));
if isempty(agl10dbz)
    cloudParams.max10dbzAgl=nan;
else
    sortedAgl10=sort(agl10dbz,'descend');
    sortedAgl10(isnan(sortedAgl10))=[];
    percIndAgl10=ceil(percWanted*length(sortedAgl10));
    cloudParams.max10dbzAgl=sortedAgl10(percIndAgl10);
end

% Cloud thickness
cloudThick=abs(maxAglAllOrig-minAglAllOrig);

cloudParams.meanThickness=mean(cloudThick,'omitnan');

sortedThick=sort(cloudThick,'descend');
sortedThick(isnan(sortedThick))=[];
percIndThick=ceil(percWanted*length(sortedThick));
cloudParams.maxThickness=sortedThick(percIndThick);

cloudThickTest=maxAglAll-minAglAll;
cloudThickTest(isnan(cloudThickTest))=[];

if length(cloudThickTest)<=length(cloudThick)/2
    cloudParams.meanThickness=nan;
end

if length(cloudThickTest)<=length(cloudThick)/2 & cloudParams.maxThickness<=10
    cloudParams.maxThickness=nan;
end

% Inhomogeneity
maxReflLin=10.^(maxRefl./10);
cloudParams.inhomo=std(maxReflLin,'omitnan')/mean(maxReflLin,'omitnan');

end