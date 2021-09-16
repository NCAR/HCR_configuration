function velDeAlias=dealiasArea(velFolded,elev)
% Unfold velocities

%foldThresh=7.8311;

velDeAlias=nan(size(velFolded));

velFolded(:,elev<0)=-velFolded(:,elev<0);

% Split in half at zero: updrafts are 1, downdrafts are 0
velHalf=nan(size(velFolded));
velHalf(velFolded>=0)=1;
velHalf(velFolded<0)=0;

% Fill nans with zeros
velHalfNoNan=velHalf;
velHalfNoNan(isnan(velFolded))=0;

% Go through updraft areas and check if they are folded
posAreas=bwconncomp(velHalfNoNan);

for ii=1:posAreas.NumObjects
    thisMask=zeros(size(velFolded));
    thisMask(posAreas.PixelIdxList{ii})=1;
    
    % Find inner boundaries
    [B,L]= bwboundaries(thisMask);
    % Outer boundaries
    edgeOut=edge(thisMask,'Sobel');
end

% velGray=mat2gray(velFolded,[-foldThresh foldThresh]);
% 
% velAdjust=imadjust(velGray,[0.4 0.6]);
% velAdjust(isnan(velFolded))=nan;
% 
% [edgeOut thresh]=edge(velAdjust,'Sobel',0.4);
% 
% velHigh=velFolded>6.5;
end