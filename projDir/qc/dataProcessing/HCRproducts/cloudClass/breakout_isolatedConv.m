function cloudPuzzleOut=breakout_isolatedConv(data)
cloudPuzzleOut=nan(size(data.DBZ));

% Make sure same pixels are non-nan
cloudPuzzleTemp=data.cloudPuzzle;
cloudPuzzleTemp(isnan(data.ECHO_TYPE_2D))=nan;

% Set conv only and elevated to mixed
echoTypeTemp=data.ECHO_TYPE_2D;
echoTypeTemp(echoTypeTemp==30 | echoTypeTemp==32)=25;

uClouds=unique(cloudPuzzleTemp(~isnan(cloudPuzzleTemp)));

newInd=1;

% Loop through clouds
for ii=1:length(uClouds)
    cloudInds=find(cloudPuzzleTemp==uClouds(ii));
    
    % Process only large regions
    if length(cloudInds)<100000
        cloudPuzzleOut(cloudInds)=newInd;
        newInd=newInd+1;
        continue
    end

    convCheck=echoTypeTemp(cloudInds);

    % Find mature convective
    % Calculate convective fraction
    convPix=sum(sum(convCheck>32));
    stratPix=sum(sum(convCheck<20));

    convFrac=convPix/(stratPix+convPix);

    % Do not process young or stratiform
    if convFrac>0.5
        cloudPuzzleOut(cloudInds)=newInd;
        newInd=newInd+1;
        continue
    end

    % Find convective areas
    % Create small mat of cloud only
    [clR,clC]=ind2sub(size(data.DBZ),cloudInds);

    typeMapBig=nan(size(data.DBZ));
    typeMapBig(cloudInds)=echoTypeTemp(cloudInds);

    typeMap=typeMapBig(min(clR):max(clR),min(clC):max(clC));

    typeMask=typeMap>25;

    % Loop through convective regions
    convRegs=bwconncomp(typeMask);

    returnInds=nan(size(typeMap));
    returnInds(~isnan(typeMap))=newInd;
    newInd=newInd+1;

    for jj=1:convRegs.NumObjects

        thisReg=zeros(size(typeMask));
        thisReg(convRegs.PixelIdxList{jj})=1;

        thisRegPadded=padarray(thisReg,[10 10],0,'both');
        typeMapPadded=padarray(typeMap,[10 10],nan,'both');

        growReg=imdilate(thisRegPadded,strel('disk',30));
        boundPix=bwperim(growReg);

        typeBound=typeMapPadded(boundPix==1);

        nanFrac=sum(isnan(typeBound))/length(typeBound);

        if nanFrac>0.7 & sum(~isnan(typeBound))<1500
            growReg=growReg(11:end-10,11:end-10);
            returnInds(growReg==1 & ~isnan(typeMap))=newInd;
            newInd=newInd+1;
        end
    end

    % Check if identified regions are still contiguous
    minObjects=min(returnInds(:),[],'omitnan');
    numObjects=max(returnInds(:),[],'omitnan');

    for cc=minObjects:numObjects
        labelMask=returnInds==cc;
        labelRegs=bwconncomp(labelMask);
        if labelRegs.NumObjects>1
            regSizes=regionprops('table',labelMask,'Area');
            sizeFrac=regSizes.Area./max(regSizes.Area);
            for dd=1:length(sizeFrac)
                if sizeFrac(dd)<1
                    thisReg=zeros(size(labelMask));
                    thisReg(labelRegs.PixelIdxList{dd})=1;
                    largerReg=imdilate(thisReg,strel('disk',2));
                    uLabSmall=unique(returnInds(thisReg==1));
                    borderVals=returnInds(largerReg==1);
                    borderVals(isnan(borderVals))=[];
                    uLab=unique(borderVals);
                    countUnique=nan(size(uLab));
                    for hh=1:length(uLab)
                        countUnique(hh)=length(find(borderVals==uLab(hh)));
                    end
                    countUnique(uLab==uLabSmall)=[];
                    uLab(uLab==uLabSmall)=[];
                    if ~isempty(uLab)
                        returnInds(labelRegs.PixelIdxList{dd})=uLab(find(countUnique==max(countUnique)));
                    end
                end
            end

        end
        cloudPuzzleOut(cloudInds)=returnInds(~isnan(returnInds));
    end
end