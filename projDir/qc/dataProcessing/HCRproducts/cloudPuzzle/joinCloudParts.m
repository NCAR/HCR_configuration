function waterMasked = joinCloudParts(waterShed)
% Join areas that have a long border or where one area is small, etc.

% Output
waterMasked=zeros(size(waterShed));
waterMasked(waterShed>0)=1;
waterMasked=bwareaopen(waterMasked,500);

waterRidges=zeros(size(waterShed));
waterRidges(waterShed==0)=1;

% Clean it up
largerRidges=imdilate(waterRidges, strel('disk', 3));

%% In the first round, we only take care of areas with large boundaries
ridgesAll=bwconncomp(largerRidges);
ridges=ridgesAll.PixelIdxList;

% Label areas with individual numbers
forLabel=bwconncomp(waterMasked);
maskLabel=labelmatrix(forLabel);            

for kk=1:ridgesAll.NumObjects
    thisRidge=ridges{kk};
            
    % Find bordering regions
    maskThis=zeros(size(largerRidges));
    maskThis(thisRidge)=1;
    
    % Skeleton
    thisRidgeSkel=bwskel(logical(maskThis));    
    ridgePixSum=sum(sum(thisRidgeSkel));
    
    if ridgePixSum>50
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
        continue
    end
        
    surrPix=maskLabel(thisRidge);
    surrPix(surrPix==0)=[];
    unPix=unique(surrPix);    
       
    if length(unPix)==1
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
        continue
    end
    
    % Check circumfirence
    cir=[];
    
    for ll=1:length(unPix)
        % Mask surrounding areas
        maskUn=zeros(size(waterShed));
        maskUn(maskLabel==unPix(ll))=1;
                
        % Circumference
        cirAll=bwconncomp(maskUn);
        cir=[cir length(cirAll.PixelIdxList{1})];
    end
    
    if length(unPix)>2
        error('Ridge divides more than one area.');
        disp(num2str(length(unPix)));
    end
    
    maxFrac=max(ridgePixSum./cir);
    if maxFrac>0.005
        waterMasked(ridges{kk})=1;
        largerRidges(ridges{kk})=1;
    end
end

%% In the second round, we only take care of areas that are too small
ridgesAll=bwconncomp(largerRidges);
ridges=ridgesAll.PixelIdxList;

% Label areas with individual numbers
forLabel=bwconncomp(waterMasked);
maskLabel=labelmatrix(forLabel);            

for kk=1:ridgesAll.NumObjects
    thisRidge=ridges{kk};
            
    % Find bordering regions
    maskThis=zeros(size(largerRidges));
    maskThis(thisRidge)=1;
    
    % Skeleton
    thisRidgeSkel=bwskel(logical(maskThis));    
    ridgePixSum=sum(sum(thisRidgeSkel));
    
    if ridgePixSum>50
        waterMasked(ridges{kk})=1;
        continue
    end
        
    surrPix=maskLabel(thisRidge);
    surrPix(surrPix==0)=[];
    unPix=unique(surrPix);    
       
    if length(unPix)==1
        waterMasked(ridges{kk})=1;
        continue
    end
    
    % Check siae of adjacent regions
    areaSize=[];
    
    for ll=1:length(unPix)
        % Mask surrounding areas
        maskUn=zeros(size(waterShed));
        maskUn(maskLabel==unPix(ll))=1;
        areaSize=[areaSize length(find(maskUn==1))];
    end
    
    if length(unPix)>2
        error('Ridge divides more than one area.');
        disp(num2str(length(unPix)));
    end
    
    minSize=min(areaSize);
    if minSize<=10000
        waterMasked(ridges{kk})=1;
    end
end
end

