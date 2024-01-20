function [velLayer,velMatRM]=buildLayer(velMat,velLayer,midPoint,highPoint,lowPoint)

velMatRM=velMat;

for kk=midPoint:0.5:highPoint
    velDiff=velMat-kk;

    [minDiff,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

    subRow=repmat((1:size(velMat,1))',[1,size(velMat,2)]);
    subCol=repmat((1:size(velMat,2)),[size(velMat,1),1]);

    minDiffIndLin=sub2ind(size(velMat),subRow,subCol,minDiffInd);

    testLayer=velMat(minDiffIndLin);
    testLayer(minDiff>0.25)=nan;

    minDiffMat=repmat(minDiff,1,1,size(velMat,3));
    minDiffLin=find(minDiffMat<=0.25);

    baseLayerMat=repmat(velLayer,1,1,size(velMat,3));
    baseLin=find(isnan(baseLayerMat));

    rmInds=intersect(minDiffIndLin,minDiffLin);
    rmInds=intersect(rmInds,baseLin);
    velMatRM(rmInds)=nan;

    velLayer(~isnan(testLayer) & isnan(velLayer))=testLayer(~isnan(testLayer) & isnan(velLayer));
end

for kk=lowPoint:0.5:midPoint
    velDiff=velMat-kk;

    [minDiff,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

    subRow=repmat((1:size(velMat,1))',[1,size(velMat,2)]);
    subCol=repmat((1:size(velMat,2)),[size(velMat,1),1]);

    minDiffIndLin=sub2ind(size(velMat),subRow,subCol,minDiffInd);

    testLayer=velMat(minDiffIndLin);
    testLayer(minDiff>0.25)=nan;

    minDiffMat=repmat(minDiff,1,1,size(velMat,3));
    minDiffLin=find(minDiffMat<=0.25);

    baseLayerMat=repmat(velLayer,1,1,size(velMat,3));
    baseLin=find(isnan(baseLayerMat));

    rmInds=intersect(minDiffIndLin,minDiffLin);
    rmInds=intersect(rmInds,baseLin);
    velMatRM(rmInds)=nan;

    velLayer(~isnan(testLayer) & isnan(velLayer))=testLayer(~isnan(testLayer) & isnan(velLayer));
end
end