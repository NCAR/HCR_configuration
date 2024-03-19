function [velLayers,powLayers]=sortMultiVelLayers(majorVel,majorPow)
velHigh=nan(size(majorVel,1),size(majorVel,2));
velLow=nan(size(majorVel,1),size(majorVel,2));
powHigh=nan(size(majorVel,1),size(majorVel,2));
powLow=nan(size(majorVel,1),size(majorVel,2));

% Find regions with more than one layer
velNums=sum(~isnan(majorVel),3);
multiLayers=nan(size(velNums));
multiLayers(velNums>1)=1;

[~,multiDens]=fast_nd_mean(multiLayers,[11,11]);

multiMask=multiDens>30;
multiMask=imclose(multiMask,strel('disk',3));
multiMask=imopen(multiMask,strel('disk',1));
multiMask=bwareaopen(multiMask,51);

mask3D=repmat(multiMask,1,1,size(majorVel,3));
majorVel(mask3D==0)=nan;
majorPow(mask3D==0)=nan;

% Find regions that are too large and break them up
multiMaskOrig=multiMask;
multiRegs=bwconncomp(multiMask);

bigInds=find((cellfun(@length, multiRegs.PixelIdxList))>1000); % Default 1000
indsComp=0;
regsPrev=-999;

while ~isempty(bigInds) & indsComp==0
    for kk=1:length(bigInds)
        pixLarge=multiRegs.PixelIdxList{bigInds(kk)};
        justOne=zeros(size(multiMask));
        justOne(pixLarge)=1;

        D=-bwdist(~justOne);
        mask1=imextendedmin(D,2);
        D2=imimposemin(D,mask1);
        Ld2=watershed(D2);
        justOne(Ld2==0)=0;
        multiMask(pixLarge)=justOne(pixLarge);
    end
    multiRegs=bwconncomp(multiMask);
    bigInds=find((cellfun(@length, multiRegs.PixelIdxList))>1000);
    if isequal(bigInds,regsPrev)
        indsComp=1;
    end
    regsPrev=bigInds;
end

L=labelmatrix(multiRegs);
L=double(L);
L(multiMaskOrig==0)=-99;
L(L==0)=nan;

L=fillmissing2(L,'nearest','MissingLocations',isnan(L));

% Process region by region
for ll=1:max(L(:))
    pixLin=find(L==ll);

    % Check how big and separate
    [r,c]=ind2sub(size(multiMask),pixLin);
    minR=min(r);
    maxR=max(r);
    minC=min(c);
    maxC=max(c);

    justOne=zeros(size(multiMask));
    justOne(pixLin)=1;

    justOne=justOne(minR:maxR,minC:maxC);
   
    [velHigh(pixLin),velLow(pixLin),powHigh(pixLin),powLow(pixLin)]=processMultiVelRegs(majorVel(minR:maxR,minC:maxC,:),majorPow(minR:maxR,minC:maxC,:), ...
        multiLayers(minR:maxR,minC:maxC),justOne);
end

velLayers=cat(3,velLow,velHigh);
powLayers=cat(3,powLow,powHigh);

end