function maskLast = thresholdMask(inData,startThresh,minArea)
% Find optimal reflectivity threshold for clouds
minRefl=min(min(inData));

maskThresh=zeros(size(inData));
maskThresh(inData>startThresh)=1;

maskFilled= imfill(maskThresh,'holes');
maskBig=bwareaopen(maskFilled,500);

maskConn=bwconncomp(maskBig);
maskNum=maskConn.NumObjects;

outNum=maskNum;
newThresh=startThresh-1;

maskGoodL=maskBig;

if newThresh<=minRefl
    maskLast=maskGoodL;
end

% Enlarge thresholding mask as long as it has the same number of
% objects
while outNum==maskNum & newThresh>minRefl
    maskLast=maskGoodL;
    
    % Make bigger mask
    maskThreshL=zeros(size(inData));
    maskThreshL(inData>newThresh)=1;
    
    maskFilledL= imfill(maskThreshL,'holes');
    maskBigL=bwareaopen(maskFilledL,minArea);
    
    maskConnL=bwconncomp(maskBigL);
    maskNumL=maskConnL.NumObjects;
    
    maskGoodL=zeros(size(maskBig));
    
    % Find new areas that overlap old areas
    if maskNumL>1
        for ii=1:maskNumL
            maskTest=maskGoodL;
            compMask=zeros(size(maskBig));
            area=maskConnL.PixelIdxList{ii};
            compMask(area)=1;
            maskAdded=compMask+maskBig;
            if max(max(maskAdded))>1
                maskGoodL(area)=1;
            end
        end
    end
    
    maskConnEnd=bwconncomp(maskGoodL);
    outNum=maskConnEnd.NumObjects;
    
    newThresh=newThresh-1;
end
end

