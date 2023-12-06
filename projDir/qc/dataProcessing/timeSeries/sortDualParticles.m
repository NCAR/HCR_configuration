function momentsVelDual=sortDualParticles(momentsVelDualRaw,momentsTime)
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
regsHigh=bwareaopen(regsHigh,25);

for ii=1:size(velHigh,3)
    thisHigh=velHigh(:,:,ii);
    thisHigh(regsHigh==0)=nan;
    velHigh(:,:,ii)=thisHigh;
end

velHigh=sort(velHigh,3);
findEmptyHigh=squeeze(sum(sum(~isnan(velHigh),1),2));
velHigh(:,:,findEmptyHigh==0)=[];

regsLow=any(~isnan(velLow),3);
regsLow=bwareaopen(regsLow,25);

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

% Add base layer to high and low
velHighFilled=velHighHoles;
velHighFilled(isnan(velHighHoles))=baseLayer(isnan(velHighHoles));

velLowFilled=velLowHoles;
velLowFilled(isnan(velLowHoles))=baseLayer(isnan(velLowHoles));

momentsVelDual(:,:,1)=baseLayer;
momentsVelDual(:,:,2)=velHighFilled;
momentsVelDual(:,:,3)=velLowFilled;
momentsVelDual(:,:,4)=velHighHoles;
momentsVelDual(:,:,5)=velLowHoles;
end