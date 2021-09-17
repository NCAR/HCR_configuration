function velDeAlias=dealiasAreaNeg(velFolded,nyq)
% Unfold velocities
velFolded=-velFolded;
velDeAlias=velFolded;

% Split in half at zero: updrafts are 1, downdrafts are 0
velFoldedTest=velFolded;
velFoldedTest(velFolded>-2 & velFolded<2)=0;
velHalf=nan(size(velFolded));
velHalf(velFoldedTest>0)=1;
velHalf(velFoldedTest<=0)=0;

% Fill nans with zeros
velHalfNoNan=velHalf;
velHalfNoNan(isnan(velFolded))=0;

% Go through updraft areas and check if they are folded
posAreas=bwconncomp(velHalfNoNan);

for ii=1:posAreas.NumObjects
    
    thisMask=zeros(size(velFolded));
    thisMask(posAreas.PixelIdxList{ii})=1;
    
    % Find boundaries
    [B,L]= bwboundaries(thisMask);
    
    % Enlarge first boundary
    boundMask=zeros(size(velFolded));
    boundMask(sub2ind(size(boundMask),B{1}(:,1),B{1}(:,2)))=1;
    
    largeBound=imdilate(boundMask,strel('disk',1));
    
    % Velocity of large boundary
    velLbound=velDeAlias(largeBound==1);
    
    % Small areas: check value of minimum and maximum difference
    if sum(~isnan(velLbound))<=5
        maxDiff=max(velLbound)-min(velLbound);
        velThis=velDeAlias(thisMask==1);
        % Dealias
        if maxDiff>14 & max(velThis)>7
            velDeAlias(thisMask==1)=velDeAlias(thisMask==1)-(2*nyq);
        end
    else
        
        % Large areas: if standard deviation is large, it is a folded area
        stdVel=std(velLbound,'omitnan');
        
        if stdVel>5
            velDeAlias(thisMask==1)=velDeAlias(thisMask==1)-(2*nyq);
        end
    end
end

velDeAlias=-velDeAlias;
end