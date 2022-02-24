function [velDeAlias doNeg]=dealiasAreaLarge(velFolded,nyq)
% Unfold velocities
velDeAlias=velFolded;
% doNeg=0;
% 
% % Split in half at zero: updrafts are 1, downdrafts are 0
% velFoldedTest=velFolded;
% velFoldedTest(velFolded>-1 & velFolded<1)=0;
% velHalf=nan(size(velFolded));
% velHalf(velFoldedTest>0)=1;
% velHalf(velFoldedTest<=0)=0;
% 
% % Fill nans with zeros
% velHalfNoNan=velHalf;
% velHalfNoNan(isnan(velFolded))=0;

% Check if any folding occurs
diffHor=diff(velFolded,1,2);
diffHor=diffHor(1:end-1,:);
diffVer=diff(velFolded,1,1);
diffVer=diffVer(:,1:end-1);

diffBoth=zeros(size(diffVer));

diffBoth(abs(diffHor)>(nyq+nyq/2))=1;
diffBoth(abs(diffVer)>(nyq+nyq/2))=1;

diffLines=bwconncomp(diffBoth);
sizeLines=sort(cellfun(@length,diffLines.PixelIdxList));

if max(sizeLines)<=5
    return
end

% Double positive folding
doublePosMask=zeros(size(velFolded));
doublePosMask(velFolded>nyq*2-1)=1;

doublePosMask=imfill(doublePosMask,'holes');

velDeAlias=unfold(doublePosMask,velDeAlias,nyq,@minus);

% Other positives
oncePosMask=zeros(size(velFolded));
oncePosMask(velDeAlias>nyq-1)=1;

oncePosMask=imfill(oncePosMask,'holes');

velDeAlias=unfold(oncePosMask,velDeAlias,nyq,@minus);

% Double negative folding
doubleNegMask=zeros(size(velFolded));
doubleNegMask(velDeAlias<-(nyq*2-1))=1;

doubleNegMask=imfill(doubleNegMask,'holes');

velDeAlias=unfold(doubleNegMask,velDeAlias,nyq,@plus);

% Other Negatives
onceNegMask=zeros(size(velFolded));
onceNegMask(velDeAlias<-(nyq-1))=1;

onceNegMask=imfill(onceNegMask,'holes');

velDeAlias=unfold(onceNegMask,velDeAlias,nyq,@plus);

end