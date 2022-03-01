function [velDeAlias,diffOut]=unfoldRegions(diffBoth,velDeAlias,nyq)
% Unfold areas that are folded
diffOut=zeros(size(diffBoth));

% Create zero regions
zeroMask=isnan(velDeAlias);
zeroMask(velDeAlias<1)=1;

zeroMask2=imdilate(zeroMask,strel('disk',3));
zeroMask2=bwskel(zeroMask2,'MinBranchLength',20);

zeroMask(zeroMask2==1)=1;

areaMask=abs(zeroMask-1);
areaMask=imfill(areaMask,'holes');

% Prepare outer boundary mask
noVelMask=isnan(velDeAlias);
noVelMask=imdilate(noVelMask,strel('disk',1));

% Go through areas and check if they are folded
singAreas=bwconncomp(areaMask,4);

for ii=1:singAreas.NumObjects

    % Check for high values
    velThis=velDeAlias(singAreas.PixelIdxList{ii});
    isFolded=length(find(velThis>nyq+nyq/2));

    if isFolded<5
        continue
    end

    thisMask=zeros(size(zeroMask));
    thisMask(singAreas.PixelIdxList{ii})=1;

    % Find boundaries
    [B,L]= bwboundaries(thisMask);

    % Enlarge first boundary
    boundMask=zeros(size(zeroMask));
    boundMask(sub2ind(size(boundMask),B{1}(:,1),B{1}(:,2)))=1;

    boundMask(noVelMask==1)=0;

    sumBound=sum(sum(boundMask));

    largeBound=imdilate(boundMask,strel('disk',1));

    diffBound=diffBoth(largeBound==1);

    sumFold=sum(diffBound);

    foldFrac=sumFold/sumBound;

    if foldFrac>0.1 % This may be too large
        velDeAlias(thisMask==1)=velDeAlias(thisMask==1)-(2*nyq);
        diffOut(largeBound==1)=1;
    end

%     % Velocity of large boundary
%     velLbound=velDeAlias(largeBound==1);
% 
%     dealias=0;
% 
%     % Small areas: check value of minimum and maximum difference
%     if sum(~isnan(velLbound))<=5
%         maxDiff=max(velLbound)-min(velLbound);
%         velThis=velDeAlias(thisMask==1);
%         % Dealias
%         if maxDiff>14 & max(velThis)>7
%             velDeAlias(thisMask==1)=plusMinus(velDeAlias(thisMask==1),(2*nyq));
%             dealias=1;
% 
%             % Check if it is working
%             velLboundC=velDeAlias(largeBound==1);
% 
%             maxDiffC=max(velLboundC)-min(velLboundC);
%             velThisC=velDeAlias(thisMask==1);
%             % Dealias back
%             if maxDiffC>5 & max(velThisC)>10
%                 velDeAlias(thisMask==1)=velDeAlias(thisMask==1);
%                 dealias=0;
%                 doNeg=1;
%             end
%         end
%     else
% 
%         % Large areas: if standard deviation is large, it is a folded area
%         stdVel=std(velLbound,'omitnan');
% 
%         if stdVel>5
%             velDeAlias(thisMask==1)=plusMinus(velDeAlias(thisMask==1),(2*nyq));
%             dealias=1;
% 
%             % Check if it is working
%             velLboundC=velDeAlias(largeBound==1);
%             stdVelC=std(velLboundC,'omitnan');
% 
%             if stdVelC>3
%                 velDeAlias(thisMask==1)=velDeAlias(thisMask==1);
%                 dealias=0;
%                 doNeg=1;
%             end
%         end
%     end

end
end