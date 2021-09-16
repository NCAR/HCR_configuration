function velDeAlias=dealiasArea(velFolded,elev)
% Unfold velocities

velFolded(:,elev>0)=-velFolded(:,elev>0);
velDeAlias=velFolded;

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
    
    if rem(ii,100)==0
        disp(['Checking area ',num2str(ii),' of ',num2str(posAreas.NumObjects),' ...']);
    end
    
    thisMask=zeros(size(velFolded));
    thisMask(posAreas.PixelIdxList{ii})=1;
    
    % Find boundaries
    [B,L]= bwboundaries(thisMask);
        
    % Enlarge first boundary
    boundMask=zeros(size(velFolded));
    boundMask(sub2ind(size(boundMask),B{1}(:,1),B{1}(:,2)))=1;
   
    largeBound=imdilate(boundMask,strel('disk',1));
    
    % Velocity of large boundary
    velLbound=velFolded(largeBound==1);
    
    dealias=0;
    
    % Small areas: check value of minimum and maximum difference
    if sum(~isnan(velLbound))<=5
%         maxDiff=max(velLbound)-min(velLbound);
%         velThis=velFolded(thisMask==1);
%         % Dealias
%         if maxDiff>14 & max(velThis)>7
%             velDeAlias(thisMask==1)=velFolded(thisMask==1)-16;
%             dealias=1;
%         end
        continue
    end
    
    % Large areas: if standard deviation is large or border has lots of large values, it is a folded area
    stdVel=std(velLbound,'omitnan');
    largeBound=sum(abs(velFolded(thisMask==1))>12);
        
    if stdVel>4 & largeBound>5
        velDeAlias(thisMask==1)=velFolded(thisMask==1)-16;
        dealias=1;
    end
    
%     % Check for holes
%     if length(B)>1
%         if dealias==0
%             for jj=2:length(B)
%                 % Enlarge boundary
%                 boundMask=zeros(size(velFolded));
%                 boundMask(sub2ind(size(boundMask),B{jj}(:,1),B{jj}(:,2)))=1;
%                 
%                 largeBound=imdilate(boundMask,strel('disk',1));
%                 
%                 % Velocity of large boundary
%                 velLbound=velFolded(largeBound==1);
%                                 
%                 % Small areas: check value of minimum and maximum difference
%                 if sum(~isnan(velLbound))<=5
% %                     maxDiff=max(velLbound)-min(velLbound);
% %                     velThis=velFolded(L==jj);
% %                     % Dealias
% %                     if maxDiff>14 & max(velThis)<-7
% %                         velDeAlias(L==jj)=velFolded(L==jj)+16;
% %                     end
%                     continue
%                 end
%                 
%                 % Large areas: if standard deviation is large, it is a folded area
%                 stdVel=std(velLbound,'omitnan');
%                 
%                 if stdVel>5
%                     velDeAlias(L==jj)=velFolded(L==jj)+16;
%                 end
%             end            
%         end
%     end
end

velDeAlias(:,elev>0)=-velDeAlias(:,elev>0);

end