function [velLayer,powLayer,velMatRM,powMatRM]=buildLayer(velMat,velLayer,powMat,powLayer,midPoint,highPoint,lowPoint)

velMatRM=velMat;
powMatRM=powMat;

stepsUp=midPoint:0.5:highPoint+0.5;
stepsDown=lowPoint-0.5:0.5:midPoint;
stepsDown=stepsDown-(stepsUp(1)-stepsDown(end));
stepsDown=fliplr(stepsDown);
maxLength=max([length(stepsUp),length(stepsDown)]);

stepsUN=nan(1,maxLength);
stepsDN=stepsUN;
stepsUN(1:length(stepsUp))=stepsUp;
stepsDN(1:length(stepsDown))=stepsDown;
steps=[stepsUN;stepsDN];
steps=steps(:)';
steps(isnan(steps))=[];

for ll=1:length(steps)
    kk=steps(ll);
    
    velDiff=velMat-kk;

    [minDiff,minDiffInd]=min(abs(velDiff),[],3,'omitmissing');

    subRow=repmat((1:size(velMat,1))',[1,size(velMat,2)]);
    subCol=repmat((1:size(velMat,2)),[size(velMat,1),1]);

    minDiffIndLin=sub2ind(size(velMat),subRow,subCol,minDiffInd);

    testLayer=velMat(minDiffIndLin);
    testLayer(minDiff>0.25)=nan;

    testLayerP=powMat(minDiffIndLin);
    testLayerP(minDiff>0.25)=nan;

    minDiffMat=repmat(minDiff,1,1,size(velMat,3));
    minDiffLin=find(minDiffMat<=0.25);

    baseLayerMat=repmat(velLayer,1,1,size(velMat,3));
    baseLin=find(isnan(baseLayerMat));

    rmInds=intersect(minDiffIndLin,minDiffLin);
    rmInds=intersect(rmInds,baseLin);
    velMatRM(rmInds)=nan;
    powMatRM(rmInds)=nan;

    powLayer(~isnan(testLayer) & isnan(velLayer))=testLayerP(~isnan(testLayer) & isnan(velLayer));
    velLayer(~isnan(testLayer) & isnan(velLayer))=testLayer(~isnan(testLayer) & isnan(velLayer));
end
end