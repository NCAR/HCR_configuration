function [noiseThresh,meanNoiseOut,R2]=findNoiseThreshMat(powInAll,avNum)
% Find noise threshold and mean noise following
% Hildebrand and Sekhon, 1974 https://doi.org/10.1175/1520-0450(1974)013%3C0808:ODOTNL%3E2.0.CO;2

%% Third version

noiseThresh=nan(size(powInAll,1),1);
meanNoiseOut=nan(size(powInAll,1),1);
R2=nan(size(powInAll,1),1);

goodInds=find(any(~isnan(powInAll),2));
powIn=powInAll(goodInds,:);

powLin2=10.^(powIn./10);

testThresh=sort(powLin2,2,'descend');
testThresh=repmat(testThresh,1,1,size(testThresh,2));

testPow=permute(testThresh,[1,3,2]);

testPow(testPow>testThresh)=nan;
sampleNum2=squeeze(sum(~isnan(testPow),3));
meanNoiseMat=squeeze(sum(testPow,3,'omitmissing'))./sampleNum2;
Qmat=squeeze(sum(testPow.^2./sampleNum2,3,'omitmissing'))-meanNoiseMat.^2;
R2mat=meanNoiseMat.^2./(Qmat.*avNum);

R2mask=R2mat>1;
[~,R2indCol]=max(R2mask,[],2);
R2indCol=R2indCol-1;
R2indRow=1:length(goodInds);
notFound=goodInds(find(R2indCol==0));
R2indCol(R2indCol==0)=1;
R2ind=sub2ind(size(meanNoiseMat),R2indRow',R2indCol);
meanNoiseMax2=meanNoiseMat(R2ind);
testThreshTest=squeeze(testThresh(:,:,1));
noiseThreshMax2=testThreshTest(R2ind);

noiseThresh(goodInds)=real(10.*log10(noiseThreshMax2));
meanNoiseOut(goodInds)=10.*log10(meanNoiseMax2);
R2(goodInds)=R2mat(R2ind);

if ~isempty(notFound)
    noiseThresh(notFound)=nan;
    meanNoiseOut(notFound)=nan;
    R2(notFound)=nan;
end

end