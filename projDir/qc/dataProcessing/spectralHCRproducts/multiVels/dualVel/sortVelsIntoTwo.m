function [momentsVelDual,momentsDbzDual]=sortVelsIntoTwo(momentsVelDualRaw,momentsPowDualRaw,momentsTime)
rmThresh=5; % Default 25

% Correct velocity for aircraft motion
xCorr=sind(momentsTime.azimuth_vc).*cosd(momentsTime.elevation).*momentsTime.eastward_velocity;
yCorr=cosd(momentsTime.azimuth_vc).*cosd(momentsTime.elevation).*momentsTime.northward_velocity;
zCorr=sind(momentsTime.elevation).*momentsTime.vertical_velocity;
momentsVelDualC=momentsVelDualRaw+xCorr+yCorr+zCorr;

% Sort vel dual into two fields
momentsVelDual=nan(size(momentsVelDualC,1),size(momentsVelDualC,2),7);
momentsDbzDual=nan(size(momentsVelDualC,1),size(momentsVelDualC,2),7);

velDiff=momentsVelDualC-momentsTime.vel;

% Fill in the base layer, i.e. the layer that is closest to the time domain
% vel
[~,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

subRow=repmat((1:size(momentsVelDualC,1))',[1,size(momentsVelDualC,2)]);
subCol=repmat((1:size(momentsVelDualC,2)),[size(momentsVelDualC,1),1]);

minDiffIndLin=sub2ind(size(momentsVelDualC),subRow,subCol,minDiffInd);

baseLayer=momentsVelDualC(minDiffIndLin);
baseLayerWork=baseLayer;

baseLayerDbz=momentsPowDualRaw(minDiffIndLin);
baseLayerDbzWork=baseLayerDbz;

% Handle the other layers
otherLayers=momentsVelDualC;
otherLayers(minDiffIndLin)=nan;

otherLayersDbz=momentsPowDualRaw;
otherLayersDbz(minDiffIndLin)=nan;

% Separate into higher and lower than from time domain
velHigh=otherLayers;
velHigh(velDiff<0)=nan;

velLow=otherLayers;
velLow(velDiff>=0)=nan;

dbzHigh=otherLayersDbz;
dbzHigh(velDiff<0)=nan;

dbzLow=otherLayersDbz;
dbzLow(velDiff>=0)=nan;

% Combine to remove small regions
regsHigh=any(~isnan(velHigh),3);
regsHigh=bwareaopen(regsHigh,rmThresh);

for ii=1:size(velHigh,3)
    thisHigh=velHigh(:,:,ii);
    thisHigh(regsHigh==0)=nan;
    velHigh(:,:,ii)=thisHigh;

    thisHighDbz=dbzHigh(:,:,ii);
    thisHighDbz(regsHigh==0)=nan;
    dbzHigh(:,:,ii)=thisHighDbz;
end

[d1,d2]=meshgrid(1:size(velHigh,2),1:size(velHigh,1));
[velHigh,Ih]=sort(velHigh,3);
dbzHighNan=nan(size(velHigh));
for jj=1:size(dbzHighNan,3)
    linInds=sub2ind(size(dbzHighNan),d2,d1,Ih(:,:,jj));
    dbzHighNan(:,:,jj)=dbzHigh(linInds);
end

dbzHigh=dbzHighNan;

findEmptyHigh=squeeze(sum(sum(~isnan(velHigh),1),2));
velHigh(:,:,findEmptyHigh==0)=[];
dbzHigh(:,:,findEmptyHigh==0)=[];

regsLow=any(~isnan(velLow),3);
regsLow=bwareaopen(regsLow,rmThresh);

for ii=1:size(velLow,3)
    thisLow=velLow(:,:,ii);
    thisLow(regsLow==0)=nan;
    velLow(:,:,ii)=thisLow;

    thisLowDbz=dbzLow(:,:,ii);
    thisLowDbz(regsLow==0)=nan;
    dbzLow(:,:,ii)=thisLowDbz;
end

[velLow,Il]=sort(velLow,3);
dbzLowNan=nan(size(velLow));
for jj=1:size(dbzLowNan,3)
    linInds=sub2ind(size(dbzLowNan),d2,d1,Il(:,:,jj));
    dbzLowNan(:,:,jj)=dbzLow(linInds);
end

dbzLow=dbzLowNan;

findEmptyLow=squeeze(sum(sum(~isnan(velLow),1),2));
velLow(:,:,findEmptyLow==0)=[];
dbzLow(:,:,findEmptyLow==0)=[];

if ~isempty(velHigh)
    multiSum=sum(~isnan(velHigh),3);
    multiSum=padarray(multiSum,[1,1],0,'both');

    % Get regions with only one
    highLayer=velHigh(:,:,1);
    highLayer=padarray(highLayer,[1,1],nan,'both');
    highLayer(multiSum~=1)=nan;

    highLayerDbz=dbzHigh(:,:,1);
    highLayerDbz=padarray(highLayerDbz,[1,1],nan,'both');
    highLayerDbz(multiSum~=1)=nan;

    velHigh=padarray(velHigh,[1,1,0],nan,'both');
    dbzHigh=padarray(dbzHigh,[1,1,0],nan,'both');

    % Pixels with more
    [r,c]=find(multiSum>1);

    for ii=1:length(r)
        surrPix=highLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
        surrMed=median(surrPix(:),'omitmissing');

        multPix=squeeze(velHigh(r(ii),c(ii),:));
        multPixDbz=squeeze(dbzHigh(r(ii),c(ii),:));
        [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
        highLayer(r(ii),c(ii))=multPix(goodInd);
        highLayerDbz(r(ii),c(ii))=multPixDbz(goodInd);
    end
    velHighHoles=highLayer(2:end-1,2:end-1);
    dbzHighHoles=highLayerDbz(2:end-1,2:end-1);
end

if ~isempty(velLow)
    multiSum=sum(~isnan(velLow),3);
    multiSum=padarray(multiSum,[1,1],0,'both');

    % Get regions with only one
    lowLayer=velLow(:,:,1);
    lowLayer=padarray(lowLayer,[1,1],nan,'both');
    lowLayer(multiSum~=1)=nan;

    lowLayerDbz=dbzLow(:,:,1);
    lowLayerDbz=padarray(lowLayerDbz,[1,1],nan,'both');
    lowLayerDbz(multiSum~=1)=nan;

    velLow=padarray(velLow,[1,1,0],nan,'both');
    dbzLow=padarray(dbzLow,[1,1,0],nan,'both');

    % Pixels with more
    [r,c]=find(multiSum>1);

    for ii=1:length(r)
        surrPix=lowLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
        surrMed=median(surrPix(:),'omitmissing');

        multPix=squeeze(velLow(r(ii),c(ii),:));
        multPixDbz=squeeze(dbzLow(r(ii),c(ii),:));
        [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
        lowLayer(r(ii),c(ii))=multPix(goodInd);
        lowLayerDbz(r(ii),c(ii))=multPixDbz(goodInd);
    end
    velLowHoles=lowLayer(2:end-1,2:end-1);
    dbzLowHoles=lowLayerDbz(2:end-1,2:end-1);
end

% Remove regions where velocities are close together
% velHighHoles(abs(velHighHoles-baseLayerWork)<0)=nan;
% velLowHoles(abs(velLowHoles-baseLayerWork)<0)=nan;

% Sort regions with only two peaks
% If data in base and high, move base to low
velHighFilled=velHighHoles;
velHighFilled(isnan(velHighHoles) & ~isnan(velLowHoles))=baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles));
baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles))=nan;

dbzHighFilled=dbzHighHoles;
dbzHighFilled(isnan(velHighHoles) & ~isnan(velLowHoles))=baseLayerDbzWork(isnan(velHighHoles) & ~isnan(velLowHoles));
baseLayerDbzWork(isnan(velHighHoles) & ~isnan(velLowHoles))=nan;

% If data in base and low, move base to high
velLowFilled=velLowHoles;
velLowFilled(~isnan(velHighHoles) & isnan(velLowHoles))=baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles));
baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles))=nan;

dbzLowFilled=dbzLowHoles;
dbzLowFilled(~isnan(velHighHoles) & isnan(velLowHoles))=baseLayerDbzWork(~isnan(velHighHoles) & isnan(velLowHoles));
baseLayerDbzWork(~isnan(velHighHoles) & isnan(velLowHoles))=nan;

% Sort regions with three peaks
% Check if base layer is closer to high or to low
minDiffPeaks=3;
diffHB=abs(velHighFilled-baseLayerWork);
diffLB=abs(velLowFilled-baseLayerWork);
% Move to closer one if difference is smaller than XX m/s
highInds=find(diffHB<diffLB & diffHB<minDiffPeaks & isnan(velHighFilled));
velHighFilled(highInds)=baseLayerWork(highInds);
baseLayerWork(diffHB<diffLB & diffHB<minDiffPeaks)=nan;
dbzHighFilled(highInds)=baseLayerDbzWork(highInds);
baseLayerDbzWork(diffHB<diffLB & diffHB<minDiffPeaks)=nan;

lowInds=(diffHB>=diffLB & diffLB<minDiffPeaks & isnan(velLowFilled));
velLowFilled(lowInds)=baseLayerWork(lowInds);
baseLayerWork(diffHB>=diffLB & diffLB<minDiffPeaks)=nan;
dbzLowFilled(lowInds)=baseLayerDbzWork(lowInds);
baseLayerDbzWork(diffHB>=diffLB & diffLB<minDiffPeaks)=nan;

% Sort regions with one peak
% Interpolate over holes and see which one is closer
intHigh=fillmissing2(velHighFilled,'movmedian',25);
intHighH=fillmissing(intHigh,'linear',2,'EndValues','extrap');
intHighV=fillmissing(intHighH,'linear',1,'EndValues','extrap');
intHigh=(intHighV+intHighH)./2;
intLow=fillmissing2(velLowFilled,'movmedian',15);
intLowH=fillmissing(intLow,'linear',2,'EndValues','extrap');
intLowV=fillmissing(intLowH,'linear',1,'EndValues','extrap');
intLow=(intLowV+intLowH)./2;

diffIntHB=abs(intHigh-baseLayerWork);
diffIntLB=abs(intLow-baseLayerWork);

indsHigh2=(diffIntHB<diffIntLB & isnan(velHighFilled));
velHighFilled(indsHigh2)=baseLayerWork(indsHigh2);
baseLayerWork(indsHigh2)=nan;
dbzHighFilled(indsHigh2)=baseLayerDbzWork(indsHigh2);
baseLayerDbzWork(indsHigh2)=nan;
indsLow2=(diffIntHB>=diffIntLB & isnan(velLowFilled));
velLowFilled(indsLow2)=baseLayerWork(indsLow2);
baseLayerWork(indsLow2)=nan;
dbzLowFilled(indsLow2)=baseLayerDbzWork(indsLow2);
baseLayerDbzWork(indsLow2)=nan;

% % velHighFilled=velHighHoles;
% % velHighFilled(isnan(velHighHoles))=baseLayer(isnan(velHighHoles));
% % 
% % velLowFilled=velLowHoles;
% % velLowFilled(isnan(velLowHoles))=baseLayer(isnan(velLowHoles));

momentsVelDual(:,:,1)=baseLayer;
momentsVelDual(:,:,2)=velHighFilled;
momentsVelDual(:,:,3)=velLowFilled;
momentsVelDual(:,:,4)=velHighHoles;
momentsVelDual(:,:,5)=velLowHoles;
momentsVelDual(:,:,6)=baseLayerWork;
momentsVelDual(:,:,7)=velHighFilled-velLowFilled;

momentsDbzDual(:,:,1)=baseLayerDbz;
momentsDbzDual(:,:,2)=dbzHighFilled;
momentsDbzDual(:,:,3)=dbzLowFilled;
momentsDbzDual(:,:,4)=dbzHighHoles;
momentsDbzDual(:,:,5)=dbzLowHoles;
momentsDbzDual(:,:,6)=baseLayerDbzWork;
momentsDbzDual(:,:,7)=dbzHighFilled-dbzLowFilled;
end