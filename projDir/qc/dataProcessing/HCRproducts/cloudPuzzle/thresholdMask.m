function maskLast = thresholdMask(inData)
% Find optimal reflectivity threshold for clouds
minRefl=round(min(min(inData)));
maxRefl=round(max(max(inData)));

allReflNumFrac=max([minRefl,-40]):min([maxRefl,0]);
allReflNumFrac=cat(2,allReflNumFrac',nan(length(allReflNumFrac),2));

totArea=length(find(~isnan(inData)));

for ii=2:size(allReflNumFrac,1)
    indII=size(allReflNumFrac,1)-ii+1;
    
    maskThresh=zeros(size(inData));
    maskThresh(inData>allReflNumFrac(indII,1))=1;
    
    maskFilled= imfill(maskThresh,'holes');
    maskBig=bwareaopen(maskFilled,500);
    
    maskConn=bwconncomp(maskBig);
    maskNum=maskConn.NumObjects;
    allReflNumFrac(indII,2)=maskNum;
    allReflNumFrac(indII,3)=sum(sum(maskBig))/totArea;
end

allReflNumFrac(allReflNumFrac(:,3)<0.75,:)=[];

[~,maxNum]=max(allReflNumFrac(:,2));

realThresh=allReflNumFrac(maxNum,1);

maskThreshOut=zeros(size(inData));
maskThreshOut(inData>realThresh)=1;

maskFilledOut= imfill(maskThreshOut,'holes');
maskLast=bwareaopen(maskFilledOut,500);

end

