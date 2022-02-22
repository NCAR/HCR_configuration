function [velFinal,changeMat]=postProcessDeAlias(velIn,nyquistVel)
% Post process velocity after spectral de-aliasing

velFinal=nan(size(velIn));
changeMat=zeros(size(velIn));

inMat=~isnan(velIn);

clouds=bwconncomp(inMat);

% Loop through clouds
for ii=1:clouds.NumObjects
    cloudInds=clouds.PixelIdxList{ii};

    [crows,ccols]=ind2sub(size(velIn),cloudInds);
    
    thisCloud=nan(max(crows)-min(crows)+1,max(ccols)-min(ccols)+1);

    smallInds=sub2ind(size(thisCloud),crows-min(crows)+1,ccols-min(ccols)+1);
    thisCloud(smallInds)=velIn(cloudInds);
end

end