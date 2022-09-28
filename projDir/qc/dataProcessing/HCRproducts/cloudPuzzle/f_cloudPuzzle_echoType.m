function cloudPuzzle=f_cloudPuzzle_echoType(cloudID,data)
cloudPuzzleSmallOnly=nan(size(data.DBZ));

uClouds=unique(cloudID(cloudID>0));

newInd=1;

% Loop through clouds
for ii=1:length(uClouds)
    cloudInds=find(cloudID==uClouds(ii));

    if cloudInds<3000000 % If small, don't do the large breakup step
        cloudPuzzleSmallOnly(cloudInds)=newInd;
        newInd=newInd+1;
    else % Break up large

        [clR,clC]=ind2sub(size(cloudID),cloudInds);

        reflMapBig=nan(size(cloudID));
        reflMapBig(cloudInds)=data.DBZ(cloudInds);

        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));

        reflStrong=reflMap>-15;
        reflStrong=imfill(reflStrong,'holes');

        reflStrong=bwareaopen(reflStrong,500000);

        strongAreas=bwconncomp(reflStrong);

        if strongAreas.NumObjects==1 % If only one large area with strong reflecivities
            cloudPuzzleSmallOnly(cloudInds)=newInd;
            newInd=newInd+1;
        else % If more than one are, try to break them apart
            stopHere=1;
        end
    end

end
end