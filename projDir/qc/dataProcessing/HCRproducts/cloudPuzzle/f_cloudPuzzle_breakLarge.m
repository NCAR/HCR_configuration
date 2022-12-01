function cloudPuzzle=f_cloudPuzzle_breakLarge(cloudID,data)
cloudPuzzle=nan(size(data.DBZ));

uClouds=unique(cloudID(cloudID>0));

newInd=1;

% Loop through clouds
for ii=1:length(uClouds)
    cloudInds=find(cloudID==uClouds(ii));
    if newInd==254
        stopHere=1;
    end

    if length(cloudInds)<3000000 % If small, don't break up
        cloudPuzzle(cloudInds)=newInd;
        newInd=newInd+1;
    else % Break up large
        disp('Large area found.')

        % Create small mat of cloud only
        [clR,clC]=ind2sub(size(cloudID),cloudInds);

        reflMapBig=nan(size(cloudID));
        reflMapBig(cloudInds)=data.DBZ(cloudInds);

        reflMap=reflMapBig(min(clR):max(clR),min(clC):max(clC));

        % Threshold on reflectivity
        reflStrong=reflMap>-12;
        reflStrong=imfill(reflStrong,'holes');

        reflStrong=bwareaopen(reflStrong,500000);

        strongAreas=bwconncomp(reflStrong);

        % Test vertical extent of high refl areas
        if strongAreas.NumObjects>1
            boxAreas=regionprops('table',reflStrong,'BoundingBox');
            vertDims=boxAreas.BoundingBox(:,4);
            vertFrac=vertDims./max(vertDims);

            % Remove areas that have small vertical extent compared to
            % deepest area
            for bb=1:length(vertFrac)
                if vertFrac(bb)<0.7
                    reflStrong(strongAreas.PixelIdxList{bb})=0;
                end
            end
            strongAreas=bwconncomp(reflStrong);
        end

        % Check if there is only one area left
        if strongAreas.NumObjects==1 % If only one large area with strong reflecivities
            cloudPuzzle(cloudInds)=newInd;
            newInd=newInd+1;

            disp('No break up.')
        else % If more than one area, try to break them apart
            % Repeat thresholding with decreasing refl thresholds to cover
            % maximum area withoug joining high refl regions
            origNum=strongAreas.NumObjects;
            threshLower=-13;
            while strongAreas.NumObjects>=origNum
                reflStrongTake=reflStrong;
                reflStrong=reflMap>threshLower;
                reflStrong=imfill(reflStrong,'holes');

                reflStrong=bwareaopen(reflStrong,500000);

                strongAreas=bwconncomp(reflStrong);
                threshLower=threshLower-1;
            end
            % Repeat with 10ths of dB
            threshLower=threshLower+1;
            strongAreas.NumObjects=strongAreas.NumObjects+1;
            reflStrong=reflStrongTake;
            while strongAreas.NumObjects>=origNum
                reflStrongTake=reflStrong;
                reflStrong=reflMap>threshLower;
                reflStrong=imfill(reflStrong,'holes');

                reflStrong=bwareaopen(reflStrong,500000);

                strongAreas=bwconncomp(reflStrong);
                threshLower=threshLower-0.1;
            end

            % Calculate distance of pixels that are not in regions (with
            % too low reflectivities) to nearest region
            strongAreasTake=bwconncomp(reflStrongTake);
            labelStrong=labelmatrix(strongAreasTake);

            % bwdist calculates the distance and gives the index of the
            % nearest pixel
            [distStrong,indStrong]=bwdist(reflStrongTake);
            indStrong(isnan(reflMap))=nan;
            indStrong(distStrong==0)=nan;

            addInds=find(indStrong>0);
            indNums=indStrong(addInds);

            % Join pixels to closest region
            labelAll=labelStrong;

            for aa=1:length(addInds)
                labelAll(addInds(aa))=labelStrong(indNums(aa));
            end

            numObjects=max(labelAll(:));
            labelAll(isnan(reflMap))=nan;

            % Check if identified regions are still contiguous
            for cc=1:numObjects
                labelMask=labelAll==cc;
                labelRegs=bwconncomp(labelMask);
                if labelRegs.NumObjects>1
                    regSizes=regionprops('table',labelMask,'Area');
                    sizeFrac=regSizes.Area./max(regSizes.Area);
                    for dd=1:length(sizeFrac)
                        if sizeFrac(dd)<1
                            thisReg=zeros(size(labelMask));
                            thisReg(labelRegs.PixelIdxList{dd})=1;
                            largerReg=imdilate(thisReg,strel('disk',2));
                            uLabSmall=unique(labelAll(thisReg==1));
                            uLab=unique(labelAll(largerReg==1));
                            uLab(uLab==0)=[];
                            uLab(uLab==uLabSmall)=[];
                            if length(uLab)==1
                                labelAll(labelRegs.PixelIdxList{dd})=uLab;
                            else
                                disp('Something is wrong with labels.')
                            end
                        end
                    end

                end
            end

            cloudPuzzle(cloudInds)=labelAll(labelAll>0)+newInd-1;
            newInd=newInd+double(numObjects);
            disp('Break up successful.')
        end
    end
end
end