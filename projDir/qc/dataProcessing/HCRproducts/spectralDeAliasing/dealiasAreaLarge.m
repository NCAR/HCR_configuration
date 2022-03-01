function velDeAlias=dealiasAreaLarge(velFolded,nyq)
% Unfold velocities
velDeAlias=velFolded;

% Check if folding occurs
diffBoth=findFoldBoundaries(velFolded,(nyq+nyq/2));

if sum(sum(diffBoth))<=5
    return
end

diffOutAll=zeros(size(diffBoth));

% Positive folding
% oncePosMask=zeros(size(velFolded));
% oncePosMask(velDeAlias>nyq-nyq/2)=1;
% 
% oncePosMask=imfill(oncePosMask,'holes');

[velDeAlias,diffOut]=unfoldRegions(diffBoth,velDeAlias,nyq);

% Check if folding occurs
diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
diffOutAll=diffOutAll+diffOut;
diffBoth(diffOutAll==1)=0;

if sum(sum(diffBoth))<=5
    return
end

% % Negative folding
% onceNegMask=zeros(size(velFolded));
% onceNegMask(velDeAlias<-(nyq-nyq/2))=1;
% 
% onceNegMask=imfill(onceNegMask,'holes');
velDeAliasNeg=-velDeAlias;
[velDeAliasNeg,diffOut]=unfoldRegions(diffBoth,velDeAliasNeg,nyq);
velDeAlias=-velDeAliasNeg;

%% Double folding

% % Double positive folding
% doublePosMask=zeros(size(velFolded));
% doublePosMask(velFolded>nyq*2-1)=1;
% 
% doublePosMask=imfill(doublePosMask,'holes');
% 
% [velDeAlias,diffOut]=unfold(doublePosMask,diffBoth,velDeAlias,nyq,@minus);
% 
% 
% 
% % Check if folding occurs
% diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
% diffOutAll=diffOutAll+diffOut;
% diffBoth(diffOutAll==1)=0;
% 
% if sum(sum(diffBoth))<=5
%     return
% end

% % Double negative folding
% doubleNegMask=zeros(size(velFolded));
% doubleNegMask(velDeAlias<-(nyq*2-1))=1;
% 
% doubleNegMask=imfill(doubleNegMask,'holes');
% 
% [velDeAlias,diffOut]=unfold(doubleNegMask,diffBoth,velDeAlias,nyq,@plus);
% 
% Check if folding occurs
% diffBoth=findFoldBoundaries(velDeAlias,(nyq+nyq/2));
% diffOutAll=diffOutAll+diffOut;
% diffBoth(diffOutAll==1)=0;
% 
% if sum(sum(diffBoth))<=5
%     return
% end
end