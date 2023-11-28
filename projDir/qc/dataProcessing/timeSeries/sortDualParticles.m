function momentsVelDual=sortDualParticles(momentsVelDualRaw,momentsTime)
% Sort vel dual into two fields
momentsVelDual=nan(size(momentsVelDualRaw,1),size(momentsVelDualRaw,2),3);

velDiff=momentsVelDualRaw-momentsTime.vel;

% Fill in the base layer, i.e. the layer that is closest to the time domain
% vel
[~,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

subRow=repmat((1:size(momentsVelDualRaw,1))',[1,size(momentsVelDualRaw,2)]);
subCol=repmat((1:size(momentsVelDualRaw,2)),[size(momentsVelDualRaw,1),1]);

minDiffIndLin=sub2ind(size(momentsVelDualRaw),subRow,subCol,minDiffInd);

momentsVelDual(:,:,1)=momentsVelDualRaw(minDiffIndLin);

% Handle the other layers
otherLayers=momentsVelDualRaw;
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
    momentsVelDual(:,:,2)=highLayer(2:end-1,2:end-1);
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
    momentsVelDual(:,:,3)=lowLayer(2:end-1,2:end-1);
end
end