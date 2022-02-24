function [velDeAlias doNeg]=dealiasAreaLarge(velFolded,nyq)
% Unfold velocities
velDeAlias=velFolded;

% Check if folding occurs
diffBoth=findFoldBoundaries(velFolded,(nyq+nyq/2));

if sum(sum(diffBoth))<=5
    return
end

% Double positive folding
doublePosMask=zeros(size(velFolded));
doublePosMask(velFolded>nyq*2-1)=1;

doublePosMask=imfill(doublePosMask,'holes');

[velDeAlias,diffOut]=unfold(doublePosMask,diffBoth,velDeAlias,nyq,@minus);

diffOutAll=zeros(size(diffOut));

% Check if folding occurs
diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
diffOutAll=diffOutAll+diffOut;
diffBoth(diffOutAll==1)=0;

if sum(sum(diffBoth))<=5
    return
end

% Other positives
oncePosMask=zeros(size(velFolded));
oncePosMask(velDeAlias>nyq-1)=1;

oncePosMask=imfill(oncePosMask,'holes');

[velDeAlias,diffOut]=unfold(oncePosMask,diffBoth,velDeAlias,nyq,@minus);

% Check if folding occurs
diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
diffOutAll=diffOutAll+diffOut;
diffBoth(diffOutAll==1)=0;

if sum(sum(diffBoth))<=5
    return
end

% Double negative folding
doubleNegMask=zeros(size(velFolded));
doubleNegMask(velDeAlias<-(nyq*2-1))=1;

doubleNegMask=imfill(doubleNegMask,'holes');

[velDeAlias,diffOut]=unfold(doubleNegMask,diffBoth,velDeAlias,nyq,@plus);

% Check if folding occurs
diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
diffOutAll=diffOutAll+diffOut;
diffBoth(diffOutAll==1)=0;

if sum(sum(diffBoth))<=5
    return
end

% Other Negatives
onceNegMask=zeros(size(velFolded));
onceNegMask(velDeAlias<-(nyq-1))=1;

onceNegMask=imfill(onceNegMask,'holes');

[velDeAlias,diffOut]=unfold(onceNegMask,diffBoth,velDeAlias,nyq,@plus);

end