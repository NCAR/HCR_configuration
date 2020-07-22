function waterMasked = joinCloudParts(waterShed,BW2)
% Joind areas that hava a long border
ridgeLines=waterShed==0;
ridgeLines(BW2==0)=0;

ridgeSize=bwconncomp(ridgeLines);
ridges=ridgeSize.PixelIdxList;
waterMasked=waterShed;
waterMasked(BW2==0)=0;

for kk=1:length(ridges)
    thisRidge=ridges{kk};
    thisRidgeIm=zeros(size(BW2));
    thisRidgeIm(thisRidge)=1;
    largerRidge=imdilate(thisRidgeIm, strel('disk', 1));
    areaPix=waterMasked(largerRidge==1);
    areaPix(areaPix==0)=[];
    unPix=unique(areaPix);
    
    % Check circumfirence of adjacent regions
    cir=[];
    for ll=1:length(unPix)
        maskUn=zeros(size(BW2));
        maskUn(waterMasked==unPix(ll))=1;
        cirAll=bwconncomp(maskUn);
        cir=[cir length(cirAll.PixelIdxList{1})];
    end
    maxFrac=max(length(ridges{kk})./cir);
    if length(ridges{kk})>50 | maxFrac>0.001
        waterMasked(ridges{kk})=1;
    end
end
end

