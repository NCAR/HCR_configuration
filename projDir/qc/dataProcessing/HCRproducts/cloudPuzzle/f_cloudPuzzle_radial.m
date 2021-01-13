function cloudPuzzle=f_cloudPuzzle_radial(data)
% Create cloud puzzle

cloudPuzzleOut=nan(size(data.DBZ));

refl=data.DBZ;
refl(data.FLAG>1)=nan;

%% Handle missing and NS cal

%disp('Filling missing and NS cal ...');

refl = fillMissingNScal(refl,data);

%% Identify contiguous clouds that are too small

%disp('Identifying small clouds ...');

reflMask=zeros(size(refl));
reflMask(~isnan(refl))=1;

pixCut=5000;
CC = bwconncomp(reflMask);

cloudNumOrig=nan(size(data.DBZ));
countCloud=1;

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=pixCut
        cloudPuzzleOut(area)=0;
    else
        cloudNumOrig(area)=countCloud;
        countCloud=countCloud+1;
    end
end

%% Add extinct back in

%disp('Filling extinct echo ...');

[cloudNum,reflExt]=fillExtinct(data,cloudNumOrig,refl);

%% Split up individual clouds

%disp('Splitting clouds ...');

numMax=max(reshape(cloudNum,1,[]),[],'omitnan');

cloudCount=1;

for ii=1:numMax
    
    cloudInds=find(cloudNum==ii);
    
    if length(cloudInds)>100000
        
        reflMapBig=nan(size(cloudNum));
        reflMapBig(cloudInds)=reflExt(cloudInds);
        
        [clR,clC]=ind2sub(size(cloudNum),cloudInds);
        
        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
        
        % Zero padding
        reflPadded=cat(1,nan(10,size(reflMap,2)),reflMap,nan(10,size(reflMap,2)));
        reflPadded=cat(2,nan(size(reflPadded,1),10),reflPadded,nan(size(reflPadded,1),10));
        
        BW=zeros(size(reflPadded));
        BW(~isnan(reflPadded))=1;
        
        BW2=imerode(BW, strel('disk', 5));
        BW3 = bwareaopen(BW2,1000);
        
        % Distance of each cloud pixel from reflectivity threshold mask
        D = bwdist(BW3);
        
        % Watershed is an image segmentation method that looks for
        % ridges and valleys in an image
        waterShed = watershed(D);
        
        % Watershed usually over-segments so we join areas back together
        % that share a large border or where one area is too small
        
        waterCensored=double(waterShed);
        waterCensored(~BW)=nan;
        
        waterMasked=joinCloudParts(waterCensored);
        
        maskJoined=zeros(size(BW2));
        maskJoined(waterMasked>0)=1;
        
        % Reverser zero padding
        maskJoined=maskJoined(11:end-10,11:end-10);
        
        maskBack=zeros(size(reflMapBig));
        maskBack(min(clR):max(clR),min(clC):max(clC))=maskJoined;
        
        uniqueClouds=bwconncomp(maskBack);
        
        for kk=1:uniqueClouds.NumObjects
            areaU=uniqueClouds.PixelIdxList{kk};
            cloudPuzzleOut(areaU)=cloudCount;
            cloudCount=cloudCount+1;
        end
        
    else
        cloudPuzzleOut(cloudInds)=cloudCount;
        cloudCount=cloudCount+1;
    end
end

cloudPuzzleOut(isnan(reflExt))=nan;

%% Fill in pixels that are not in small areas (i.e. not zero) that have
% reflectivities but are nan in cloudPuzzle

%disp('Filling in final pixels ...');
allReflMask=zeros(size(reflExt));
allReflMask(~isnan(reflExt))=1;
allReflMask(cloudPuzzleOut==0)=0;

puzzleMask=zeros(size(reflExt));
puzzleMask(cloudPuzzleOut>0)=1;

[oldR oldC]=find(puzzleMask==1);
[addR addC]=find(puzzleMask==0 & allReflMask==1);
idx = knnsearch([oldR oldC], [addR addC]);
nearest_OldValue = cloudPuzzleOut(sub2ind(size(cloudPuzzleOut), oldR(idx), oldC(idx)));
cloudPuzzleAttached=cloudPuzzleOut;
cloudPuzzleAttached(sub2ind(size(cloudPuzzleOut), addR, addC))=nearest_OldValue;

% Sometimes areas get attached to wrong area
cloudPuzzle=cloudPuzzleAttached;
for jj=1:cloudCount-1
    maskNumber=zeros(size(cloudPuzzleAttached));
    maskNumber(cloudPuzzleAttached==jj)=1;
    individs=bwconncomp(maskNumber);
    if individs.NumObjects>1
        indivClouds=individs.PixelIdxList;
        for ll=1:individs.NumObjects
            if length(indivClouds{ll})<1001
                maskIndiv=zeros(size(cloudPuzzleAttached));
                maskIndiv(indivClouds{ll})=1;
                maskExp=imdilate(maskIndiv, strel('disk', 2));
                pixExp=cloudPuzzleAttached(find(maskExp==1));
                pixExp(find(pixExp==0 | pixExp==jj | isnan(pixExp)))=[];
                if ~isempty(pixExp)
                    pixU=unique(pixExp);
                    cloudPuzzle(indivClouds{ll})=pixU;
                end
            end
        end
    end
end

end