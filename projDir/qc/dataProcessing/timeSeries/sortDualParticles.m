function momentsVelDual=sortDualParticles(momentsVelDualRaw,momentsTime)
rmThresh=5; % Default 25

% Correct velocity for aircraft motion
xCorr=sind(momentsTime.azimuth_vc).*cosd(momentsTime.elevation).*momentsTime.eastward_velocity;
yCorr=cosd(momentsTime.azimuth_vc).*cosd(momentsTime.elevation).*momentsTime.northward_velocity;
zCorr=sind(momentsTime.elevation).*momentsTime.vertical_velocity;
momentsVelDualC=momentsVelDualRaw+xCorr+yCorr+zCorr;

% Sort vel dual into two fields
momentsVelDual=nan(size(momentsVelDualC,1),size(momentsVelDualC,2),5);

velDiff=momentsVelDualC-momentsTime.vel;

% Fill in the base layer, i.e. the layer that is closest to the time domain
% vel
[~,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

subRow=repmat((1:size(momentsVelDualC,1))',[1,size(momentsVelDualC,2)]);
subCol=repmat((1:size(momentsVelDualC,2)),[size(momentsVelDualC,1),1]);

minDiffIndLin=sub2ind(size(momentsVelDualC),subRow,subCol,minDiffInd);

baseLayer=momentsVelDualC(minDiffIndLin);
baseLayerWork=baseLayer;

% Handle the other layers
otherLayers=momentsVelDualC;
otherLayers(minDiffIndLin)=nan;

% Separate into higher and lower than from time domain
velHigh=otherLayers;
velHigh(velDiff<0)=nan;

velLow=otherLayers;
velLow(velDiff>=0)=nan;

% Combine to remove small regions
regsHigh=any(~isnan(velHigh),3);
regsHigh=bwareaopen(regsHigh,rmThresh);

for ii=1:size(velHigh,3)
    thisHigh=velHigh(:,:,ii);
    thisHigh(regsHigh==0)=nan;
    velHigh(:,:,ii)=thisHigh;
end

velHigh=sort(velHigh,3);
findEmptyHigh=squeeze(sum(sum(~isnan(velHigh),1),2));
velHigh(:,:,findEmptyHigh==0)=[];

regsLow=any(~isnan(velLow),3);
regsLow=bwareaopen(regsLow,rmThresh);

for ii=1:size(velLow,3)
    thisLow=velLow(:,:,ii);
    thisLow(regsLow==0)=nan;
    velLow(:,:,ii)=thisLow;
end
velLow=sort(velLow,3);
findEmptyLow=squeeze(sum(sum(~isnan(velLow),1),2));
velLow(:,:,findEmptyLow==0)=[];

if ~isempty(velHigh)
    multiSum=sum(~isnan(velHigh),3);
    multiSum=padarray(multiSum,[1,1],0,'both');

    % Get regions with only one
    highLayer=velHigh(:,:,1);
    highLayer=padarray(highLayer,[1,1],nan,'both');
    highLayer(multiSum~=1)=nan;

    velHigh=padarray(velHigh,[1,1,0],nan,'both');

    % Pixels with more
    [r,c]=find(multiSum>1);

    for ii=1:length(r)
        surrPix=highLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
        surrMed=median(surrPix(:),'omitmissing');

        multPix=squeeze(velHigh(r(ii),c(ii),:));
        [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
        highLayer(r(ii),c(ii))=multPix(goodInd);
    end
    velHighHoles=highLayer(2:end-1,2:end-1);
end

if ~isempty(velLow)
    multiSum=sum(~isnan(velLow),3);
    multiSum=padarray(multiSum,[1,1],0,'both');

    % Get regions with only one
    lowLayer=velLow(:,:,1);
    lowLayer=padarray(lowLayer,[1,1],nan,'both');
    lowLayer(multiSum~=1)=nan;

    velLow=padarray(velLow,[1,1,0],nan,'both');

    % Pixels with more
    [r,c]=find(multiSum>1);

    for ii=1:length(r)
        surrPix=lowLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
        surrMed=median(surrPix(:),'omitmissing');

        multPix=squeeze(velLow(r(ii),c(ii),:));
        [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
        lowLayer(r(ii),c(ii))=multPix(goodInd);
    end
    velLowHoles=lowLayer(2:end-1,2:end-1);
end

% Remove regions where velocities are close together
velHighHoles(abs(velHighHoles-baseLayerWork)<0)=nan;
velLowHoles(abs(velLowHoles-baseLayerWork)<0)=nan;

% Sort regions with only two peaks
velHighFilled=velHighHoles;
velHighFilled(isnan(velHighHoles) & ~isnan(velLowHoles))=baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles));
baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles))=nan;

velLowFilled=velLowHoles;
velLowFilled(~isnan(velHighHoles) & isnan(velLowHoles))=baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles));
baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles))=nan;

% Sort regions with three peaks
diffHB=abs(velHighFilled-baseLayerWork);
diffLB=abs(velLowFilled-baseLayerWork);

velHighFilled(diffHB<diffLB & diffHB<2)=baseLayerWork(diffHB<diffLB & diffHB<2);
baseLayerWork(diffHB<diffLB & diffHB<2)=nan;
velLowFilled(diffHB>=diffLB & diffLB<2)=baseLayerWork(diffHB>=diffLB & diffLB<2);
baseLayerWork(diffHB>=diffLB & diffLB<2)=nan;

% Sort regions with one peak
intHigh=fillmissing2(velHighFilled,'movmedian',15);
intHigh=fillmissing2(intHigh,'linear');
intLow=fillmissing2(velLowFilled,'movmedian',15);
intLow=fillmissing2(intLow,'linear');

diffIntHB=abs(intHigh-baseLayerWork);
diffIntLB=abs(intLow-baseLayerWork);

velHighFilled(diffIntHB<diffIntLB)=baseLayerWork(diffIntHB<diffIntLB);
baseLayerWork(diffIntHB<diffIntLB)=nan;
velLowFilled(diffIntHB>=diffIntLB)=baseLayerWork(diffIntHB>=diffIntLB);
baseLayerWork(diffIntHB>=diffIntLB)=nan;

% velHighFilled=velHighHoles;
% velHighFilled(isnan(velHighHoles))=baseLayer(isnan(velHighHoles));
% 
% velLowFilled=velLowHoles;
% velLowFilled(isnan(velLowHoles))=baseLayer(isnan(velLowHoles));

momentsVelDual(:,:,1)=baseLayerWork;
momentsVelDual(:,:,2)=velHighFilled;
momentsVelDual(:,:,3)=velLowFilled;
momentsVelDual(:,:,4)=velHighHoles;
momentsVelDual(:,:,5)=velLowHoles;
end