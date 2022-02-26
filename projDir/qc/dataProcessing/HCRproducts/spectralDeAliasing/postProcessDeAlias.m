function [velFinal,changeMat]=postProcessDeAlias(velIn,nyquistVel,elev)
% Post process velocity after spectral de-aliasing

velFinal=nan(size(velIn));
changeMat=zeros(size(velIn));

% Make all updrafts negative and updrafts negative
velIn(:,elev>0)=-velIn(:,elev>0);

inMat=~isnan(velIn);
clouds=bwconncomp(inMat);

% Sort by size
[sizeCloud,indSort]=sort(cellfun(@length,clouds.PixelIdxList));
pixListSorted=clouds.PixelIdxList(indSort);
pixListSorted=fliplr(pixListSorted);
sizeCloud=fliplr(sizeCloud);

largeInd=max(find(sizeCloud>50));

% Loop through large clouds
for ii=1:largeInd
%for ii=1:length(pixListSorted)
    cloudInds=pixListSorted{ii};

    [crows,ccols]=ind2sub(size(velIn),cloudInds);
    
    thisCloud=nan(max(crows)-min(crows)+1,max(ccols)-min(ccols)+1);

    smallInds=sub2ind(size(thisCloud),crows-min(crows)+1,ccols-min(ccols)+1);
    thisCloud(smallInds)=velIn(cloudInds);

    % Dealias in the positive direction
    velDeAliasOne=dealiasAreaLarge(thisCloud,nyquistVel);
        
    velFinal(pixListSorted{ii})=velDeAliasOne(~isnan(velDeAliasOne));
end

% Transform down pointing back to original sign
velFinal(:,elev>0)=-velFinal(:,elev>0);

end