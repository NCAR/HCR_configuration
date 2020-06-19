function cloudPuzzleOut = f_cloudPuzzle_radial(data)
% Make cloud puzzle

startThresh=-20;

cloudPuzzleOut=nan(size(data.DBZ));

%% Add extinct back in

disp('Filling extinct echo ...');

refl=fillExtinct(data);

%% Handle missing and NS cal

disp('Filling missing and NS cal ...');

refl = fillMissingNScal(refl,data);

%% Smooth with convolution

% Create mask for convolution
radius=5;
numPix=radius*2+1;
[rr cc] = meshgrid(1:numPix);
cirMask = sqrt((rr-(radius+1)).^2+(cc-(radius+1)).^2)<=radius;
cirMask=double(cirMask);

% Normalize
cirMask=cirMask./(sum(reshape(cirMask,1,[])));

% Convolution
reflConv=nanconv(refl,cirMask);

%% Identify contiguous clouds that are too small

disp('Identifying small clouds ...');

reflLarge=refl;

reflMask=zeros(size(refl));
reflMask(~isnan(refl))=1;

pixCut=5000;
CC = bwconncomp(reflMask);

cloudNum=nan(size(data.DBZ));
countCloud=1;

for ii=1:CC.NumObjects
    area=CC.PixelIdxList{ii};
    if length(area)<=pixCut
        reflLarge(area)=nan;
        cloudPuzzleOut(area)=0;
    else
        cloudNum(area)=countCloud;
        countCloud=countCloud+1;
    end
end

%% Split up individual clouds

disp('Splitting clouds ...');

numMax=max(reshape(cloudNum,1,[]),[],'omitnan');

cloudCount=1;

for ii=1:numMax
    
    cloudInds=find(cloudNum==ii);
    
    if length(cloudInds)>100000
             
        reflMapBig=nan(size(cloudNum));
        reflMapBig(cloudInds)=reflConv(cloudInds);
        
        [clR clC]=ind2sub(size(cloudNum),cloudInds);
        
        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));
        aslMap=data.asl(min(clR):max(clR),min(clC):max(clC));
        timeMap=data.time(min(clC):max(clC));
        
        % Zero padding
        reflPadded=cat(1,nan(10,size(reflMap,2)),reflMap,nan(10,size(reflMap,2)));
        reflPadded=cat(2,nan(size(reflPadded,1),10),reflPadded,nan(size(reflPadded,1),10));
                     
        % BW mask with filled holes
        BW=zeros(size(reflPadded));
        BW(~isnan(reflPadded))=1;
        
        BW2 = imfill(BW,'holes');
                 
        % Distance
        D = -bwdist(~BW2);
               
        % Thresholding
        maskBig = thresholdMask(reflPadded,startThresh,500);
                
        newMin = imimposemin(D,maskBig);
        waterShed = watershed(newMin);
        
        waterMasked=joinCloudParts(waterShed,BW2);
        
        maskJoined=zeros(size(BW2));
        maskJoined(waterMasked>0)=1;
        
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
end
