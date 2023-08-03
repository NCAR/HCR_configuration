function [velDeAlias doNeg]=dealiasAreaPos(velFolded,nyq)
% Unfold velocities
velDeAlias=velFolded;
doNeg=0;

% Split in half at zero: updrafts are 1, downdrafts are 0
velFoldedTest=velFolded;
velFoldedTest(velFolded>-1 & velFolded<1)=0;
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
    
    dealias=0;
    
    % Small areas: check value of minimum and maximum difference
    if sum(~isnan(velLbound))<=5
        maxDiff=max(velLbound)-min(velLbound);
        velThis=velDeAlias(thisMask==1);
        % Dealias
        if maxDiff>14 & max(velThis)>7
            velDeAlias(thisMask==1)=velDeAlias(thisMask==1)-(2*nyq);
            dealias=1;
            
            % Check if it is working
            velLboundC=velDeAlias(largeBound==1);
            
            maxDiffC=max(velLboundC)-min(velLboundC);
            velThisC=velDeAlias(thisMask==1);
            % Dealias back
            if maxDiffC>5 & max(velThisC)>10
                velDeAlias(thisMask==1)=velDeAlias(thisMask==1);
                dealias=0;
                doNeg=1;
            end
        end
    else
        
        % Large areas: if standard deviation is large, it is a folded area
        stdVel=std(velLbound,'omitnan');
        
        if stdVel>5
            velDeAlias(thisMask==1)=velDeAlias(thisMask==1)-(2*nyq);
            dealias=1;
            
            % Check if it is working
            velLboundC=velDeAlias(largeBound==1);
            stdVelC=std(velLboundC,'omitnan');
            
            if stdVelC>3
                velDeAlias(thisMask==1)=velDeAlias(thisMask==1);
                dealias=0;
                doNeg=1;
            end
        end
    end
    
    % Check for holes
    if length(B)>1
        doNeg=1;
        for jj=2:length(B)
            % Enlarge boundary
            boundMask=zeros(size(velFolded));
            boundMask(sub2ind(size(boundMask),B{jj}(:,1),B{jj}(:,2)))=1;
            
            largeBound=imdilate(boundMask,strel('disk',1));
            
            % Velocity of large boundary
            velLbound=velDeAlias(largeBound==1);
            
            % Small areas: check value of minimum and maximum difference
            if sum(~isnan(velLbound))<=5
                maxDiff=max(velLbound)-min(velLbound);
                velThis=velDeAlias(L==jj);
                % Dealias
                if maxDiff>14 & max(velThis)<-7
                    velDeAlias(L==jj)=velDeAlias(L==jj)+(2*nyq);
                end
                continue
            end
            
            % Large areas: if standard deviation is large, it is a folded area
            stdVel=std(velLbound,'omitnan');
            
            if stdVel>5
                velDeAlias(L==jj)=velDeAlias(L==jj)+(2*nyq);
            end
        end
    end
end

end