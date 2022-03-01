function velDeAlias=dealiasAreaLarge(velFolded,nyq)
% Unfold velocities
velDeAlias=velFolded;

diffOutAll=zeros(size(velFolded));
diffOut=zeros(size(velFolded));

%% Positive

% Find number of folding
maxVel=max(reshape(velFolded,1,[]),[],'omitnan');
numFold=floor(maxVel/nyq);

for ii=1:numFold

    % Double
    velDeAlias=velDeAlias-(numFold-ii+1)*nyq;

    % Check if folded
    diffBoth=findFoldBoundaries(velDeAlias,(nyq/2));
    diffOutAll=diffOutAll+diffOut;
    diffBoth(diffOutAll==1)=0;

    if sum(sum(diffBoth))<=5
        return
    end

    % Unfold
    [velDeAlias,diffOut]=unfoldRegions(diffBoth,velDeAlias,nyq);
    velDeAlias=velDeAlias+(numFold-ii+1)*nyq;

end

%% Negative

velDeAlias=-velDeAlias;

% Find number of folding
maxVel=max(reshape(velFolded,1,[]),[],'omitnan');
numFold=floor(maxVel/nyq);

for ii=1:numFold

    % Double
    velDeAlias=velDeAlias-(numFold-ii+1)*nyq;

    % Check if folded
    diffBoth=findFoldBoundaries(velDeAlias,(nyq/2));
    diffOutAll=diffOutAll+diffOut;
    diffBoth(diffOutAll==1)=0;

    if sum(sum(diffBoth))<=5
        return
    end

    % Unfold
    [velDeAlias,diffOut]=unfoldRegions(diffBoth,velDeAlias,nyq);
    velDeAlias=velDeAlias+(numFold-ii+1)*nyq;

end

velDeAlias=-velDeAlias;

% Positive folding
% oncePosMask=zeros(size(velFolded));
% oncePosMask(velDeAlias>nyq-nyq/2)=1;
% 
% oncePosMask=imfill(oncePosMask,'holes');

% % Negative folding
% onceNegMask=zeros(size(velFolded));
% onceNegMask(velDeAlias<-(nyq-nyq/2))=1;
% 
% onceNegMask=imfill(onceNegMask,'holes');
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