function [velLayers,powLayers]=sortVelLayers(majorVel,majorPow,minorVel,minorPow,momentsTime)

% Find base layer
baseLayer=nan(size(majorVel,1),size(majorVel,2));
majorVelRM=majorVel;

for kk=3:0.5:ceil(max(majorVel(:),[],'omitmissing'))
    velDiff=majorVel-kk;

    [minDiff,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

    subRow=repmat((1:size(majorVel,1))',[1,size(majorVel,2)]);
    subCol=repmat((1:size(majorVel,2)),[size(majorVel,1),1]);

    minDiffIndLin=sub2ind(size(majorVel),subRow,subCol,minDiffInd);

    testLayer=majorVel(minDiffIndLin);
    testLayer(minDiff>0.25)=nan;

    baseLayer(~isnan(testLayer) & isnan(baseLayer))=testLayer(~isnan(testLayer) & isnan(baseLayer));
end

for kk=floor(min(majorVel(:),[],'omitmissing')):0.5:3
    velDiff=majorVel-kk;

    [minDiff,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

    subRow=repmat((1:size(majorVel,1))',[1,size(majorVel,2)]);
    subCol=repmat((1:size(majorVel,2)),[size(majorVel,1),1]);

    minDiffIndLin=sub2ind(size(majorVel),subRow,subCol,minDiffInd);

    testLayer=majorVel(minDiffIndLin);
    testLayer(minDiff>0.25)=nan;

    baseLayer(~isnan(testLayer) & isnan(baseLayer))=testLayer(~isnan(testLayer) & isnan(baseLayer));
end

surf(baseLayer,'edgecolor','none')
view(2)
caxis([-5,15])
colormap(jet(20))

% % Sort vel dual into two fields
% %velLayers=nan(size(majorVel,1),size(majorVel,2),7);
% %powLayers=nan(size(majorVel,1),size(majorVel,2),7);
% 
% velDiff=majorVel-momentsTime.vel;
% 
% % Fill in the base layer, i.e. the layer that is closest to the time domain
% % vel
% [~,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');
% 
% subRow=repmat((1:size(majorVel,1))',[1,size(majorVel,2)]);
% subCol=repmat((1:size(majorVel,2)),[size(majorVel,1),1]);
% 
% minDiffIndLin=sub2ind(size(majorVel),subRow,subCol,minDiffInd);
% 
% baseLayer=majorVel(minDiffIndLin);
% baseLayerWork=baseLayer;
% 
% 
% 
% baseLayerDbz=majorPow(minDiffIndLin);
% baseLayerDbzWork=baseLayerDbz;
% 
% % Handle the other layers
% otherLayers=majorVel;
% otherLayers(minDiffIndLin)=nan;
% 
% otherLayersDbz=majorPow;
% otherLayersDbz(minDiffIndLin)=nan;
% 
% % Separate into higher and lower than from time domain
% velHigh=otherLayers;
% velHigh(velDiff<0)=nan;
% 
% velLow=otherLayers;
% velLow(velDiff>=0)=nan;
% 
% dbzHigh=otherLayersDbz;
% dbzHigh(velDiff<0)=nan;
% 
% dbzLow=otherLayersDbz;
% dbzLow(velDiff>=0)=nan;
% 
% % Sort and remove empty layers
% [d1,d2]=meshgrid(1:size(velHigh,2),1:size(velHigh,1));
% [velHigh,Ih]=sort(velHigh,3);
% dbzHighNan=nan(size(velHigh));
% for jj=1:size(dbzHighNan,3)
%     linInds=sub2ind(size(dbzHighNan),d2,d1,Ih(:,:,jj));
%     dbzHighNan(:,:,jj)=dbzHigh(linInds);
% end
% 
% dbzHigh=dbzHighNan;
% 
% findEmptyHigh=squeeze(sum(sum(~isnan(velHigh),1),2));
% velHigh(:,:,findEmptyHigh==0)=[];
% dbzHigh(:,:,findEmptyHigh==0)=[];
% 
% [velLow,Il]=sort(velLow,3);
% dbzLowNan=nan(size(velLow));
% for jj=1:size(dbzLowNan,3)
%     linInds=sub2ind(size(dbzLowNan),d2,d1,Il(:,:,jj));
%     dbzLowNan(:,:,jj)=dbzLow(linInds);
% end
% 
% dbzLow=dbzLowNan;
% 
% findEmptyLow=squeeze(sum(sum(~isnan(velLow),1),2));
% velLow(:,:,findEmptyLow==0)=[];
% dbzLow(:,:,findEmptyLow==0)=[];
% 
% if ~isempty(velHigh)
%     multiSum=sum(~isnan(velHigh),3);
%     multiSum=padarray(multiSum,[1,1],0,'both');
% 
%     % Get regions with only one
%     highLayer=velHigh(:,:,1);
%     highLayer=padarray(highLayer,[1,1],nan,'both');
%     highLayer(multiSum~=1)=nan;
% 
%     highLayerDbz=dbzHigh(:,:,1);
%     highLayerDbz=padarray(highLayerDbz,[1,1],nan,'both');
%     highLayerDbz(multiSum~=1)=nan;
% 
%     velHigh=padarray(velHigh,[1,1,0],nan,'both');
%     dbzHigh=padarray(dbzHigh,[1,1,0],nan,'both');
% 
%     % Pixels with more
%     [r,c]=find(multiSum>1);
% 
%     for ii=1:length(r)
%         surrPix=highLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
%         surrMed=median(surrPix(:),'omitmissing');
% 
%         multPix=squeeze(velHigh(r(ii),c(ii),:));
%         multPixDbz=squeeze(dbzHigh(r(ii),c(ii),:));
%         [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
%         highLayer(r(ii),c(ii))=multPix(goodInd);
%         highLayerDbz(r(ii),c(ii))=multPixDbz(goodInd);
%     end
%     velHighHoles=highLayer(2:end-1,2:end-1);
%     dbzHighHoles=highLayerDbz(2:end-1,2:end-1);
% end
% 
% if ~isempty(velLow)
%     multiSum=sum(~isnan(velLow),3);
%     multiSum=padarray(multiSum,[1,1],0,'both');
% 
%     % Get regions with only one
%     lowLayer=velLow(:,:,1);
%     lowLayer=padarray(lowLayer,[1,1],nan,'both');
%     lowLayer(multiSum~=1)=nan;
% 
%     lowLayerDbz=dbzLow(:,:,1);
%     lowLayerDbz=padarray(lowLayerDbz,[1,1],nan,'both');
%     lowLayerDbz(multiSum~=1)=nan;
% 
%     velLow=padarray(velLow,[1,1,0],nan,'both');
%     dbzLow=padarray(dbzLow,[1,1,0],nan,'both');
% 
%     % Pixels with more
%     [r,c]=find(multiSum>1);
% 
%     for ii=1:length(r)
%         surrPix=lowLayer(r(ii)-1:r(ii)+1,c(ii)-1:c(ii)+1);
%         surrMed=median(surrPix(:),'omitmissing');
% 
%         multPix=squeeze(velLow(r(ii),c(ii),:));
%         multPixDbz=squeeze(dbzLow(r(ii),c(ii),:));
%         [~,goodInd]=min(abs(surrMed-multPix),[],'omitmissing');
%         lowLayer(r(ii),c(ii))=multPix(goodInd);
%         lowLayerDbz(r(ii),c(ii))=multPixDbz(goodInd);
%     end
%     velLowHoles=lowLayer(2:end-1,2:end-1);
%     dbzLowHoles=lowLayerDbz(2:end-1,2:end-1);
% end
% 
% % Remove regions where velocities are close together
% % velHighHoles(abs(velHighHoles-baseLayerWork)<0)=nan;
% % velLowHoles(abs(velLowHoles-baseLayerWork)<0)=nan;
% 
% % Sort regions with only two peaks
% % If data in base and high, move base to low
% velHighFilled=velHighHoles;
% velHighFilled(isnan(velHighHoles) & ~isnan(velLowHoles))=baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles));
% baseLayerWork(isnan(velHighHoles) & ~isnan(velLowHoles))=nan;
% 
% dbzHighFilled=dbzHighHoles;
% dbzHighFilled(isnan(velHighHoles) & ~isnan(velLowHoles))=baseLayerDbzWork(isnan(velHighHoles) & ~isnan(velLowHoles));
% baseLayerDbzWork(isnan(velHighHoles) & ~isnan(velLowHoles))=nan;
% 
% % If data in base and low, move base to high
% velLowFilled=velLowHoles;
% velLowFilled(~isnan(velHighHoles) & isnan(velLowHoles))=baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles));
% baseLayerWork(~isnan(velHighHoles) & isnan(velLowHoles))=nan;
% 
% dbzLowFilled=dbzLowHoles;
% dbzLowFilled(~isnan(velHighHoles) & isnan(velLowHoles))=baseLayerDbzWork(~isnan(velHighHoles) & isnan(velLowHoles));
% baseLayerDbzWork(~isnan(velHighHoles) & isnan(velLowHoles))=nan;
% 
% % Sort regions with three peaks
% % Check if base layer is closer to high or to low
% minDiffPeaks=3;
% diffHB=abs(velHighFilled-baseLayerWork);
% diffLB=abs(velLowFilled-baseLayerWork);
% % Move to closer one if difference is smaller than XX m/s
% highInds=find(diffHB<diffLB & diffHB<minDiffPeaks & isnan(velHighFilled));
% velHighFilled(highInds)=baseLayerWork(highInds);
% baseLayerWork(diffHB<diffLB & diffHB<minDiffPeaks)=nan;
% dbzHighFilled(highInds)=baseLayerDbzWork(highInds);
% baseLayerDbzWork(diffHB<diffLB & diffHB<minDiffPeaks)=nan;
% 
% lowInds=(diffHB>=diffLB & diffLB<minDiffPeaks & isnan(velLowFilled));
% velLowFilled(lowInds)=baseLayerWork(lowInds);
% baseLayerWork(diffHB>=diffLB & diffLB<minDiffPeaks)=nan;
% dbzLowFilled(lowInds)=baseLayerDbzWork(lowInds);
% baseLayerDbzWork(diffHB>=diffLB & diffLB<minDiffPeaks)=nan;
% 
% % Sort regions with one peak
% % Interpolate over holes and see which one is closer
% intHigh=fillmissing2(velHighFilled,'movmedian',25);
% intHighH=fillmissing(intHigh,'linear',2,'EndValues','extrap');
% intHighV=fillmissing(intHighH,'linear',1,'EndValues','extrap');
% intHigh=(intHighV+intHighH)./2;
% intLow=fillmissing2(velLowFilled,'movmedian',15);
% intLowH=fillmissing(intLow,'linear',2,'EndValues','extrap');
% intLowV=fillmissing(intLowH,'linear',1,'EndValues','extrap');
% intLow=(intLowV+intLowH)./2;
% 
% diffIntHB=abs(intHigh-baseLayerWork);
% diffIntLB=abs(intLow-baseLayerWork);
% 
% indsHigh2=(diffIntHB<diffIntLB & isnan(velHighFilled));
% velHighFilled(indsHigh2)=baseLayerWork(indsHigh2);
% baseLayerWork(indsHigh2)=nan;
% dbzHighFilled(indsHigh2)=baseLayerDbzWork(indsHigh2);
% baseLayerDbzWork(indsHigh2)=nan;
% indsLow2=(diffIntHB>=diffIntLB & isnan(velLowFilled));
% velLowFilled(indsLow2)=baseLayerWork(indsLow2);
% baseLayerWork(indsLow2)=nan;
% dbzLowFilled(indsLow2)=baseLayerDbzWork(indsLow2);
% baseLayerDbzWork(indsLow2)=nan;
% 
% % % velHighFilled=velHighHoles;
% % % velHighFilled(isnan(velHighHoles))=baseLayer(isnan(velHighHoles));
% % % 
% % % velLowFilled=velLowHoles;
% % % velLowFilled(isnan(velLowHoles))=baseLayer(isnan(velLowHoles));
% 
% velLayers(:,:,1)=baseLayer;
% velLayers(:,:,2)=velHighFilled;
% velLayers(:,:,3)=velLowFilled;
% velLayers(:,:,4)=velHighHoles;
% velLayers(:,:,5)=velLowHoles;
% velLayers(:,:,6)=baseLayerWork;
% velLayers(:,:,7)=velHighFilled-velLowFilled;
% 
% powLayers(:,:,1)=baseLayerDbz;
% powLayers(:,:,2)=dbzHighFilled;
% powLayers(:,:,3)=dbzLowFilled;
% powLayers(:,:,4)=dbzHighHoles;
% powLayers(:,:,5)=dbzLowHoles;
% powLayers(:,:,6)=baseLayerDbzWork;
% powLayers(:,:,7)=dbzHighFilled-dbzLowFilled;
end