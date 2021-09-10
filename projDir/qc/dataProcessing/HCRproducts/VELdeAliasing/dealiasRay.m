function velUnfold=unfoldVel(velFolded,flag,elev)
% Unfold velocities

foldThresh=7.8311;

velUnfold=nan(size(velFolded));

velFolded(flag~=1)=nan;
velFolded(:,elev<0)=-velFolded(:,elev<0);

velDiff=diff(velFolded,1,1);

goodCols=find(sum(~isnan(velFolded),1)>0);

for ii=1:length(goodCols)
    diffVec=velDiff(:,goodCols(ii));
    if abs(max(diffVec))>10
        diffVec(abs(diffVec)<=10)=nan;
        velVec=velFolded(:,goodCols(ii));
        velVecOut=velVec;
        
        % Get start end end indices
        startInds=find(diffVec>0);
        endIndsTemp=find(diffVec<0);
        if isempty(startInds) | (~isempty(endIndsTemp) & startInds(1)>endIndsTemp(1))
            startInds=[1,startInds];
        end
        if isempty(endIndsTemp)
            endIndsTemp=size(velFolded,1);
        end
        
        endInds=[];
        nextEndInd=1;
        
        if length(startInds)==1;
            endInds=endIndsTemp;
        end
        for jj=1:length(startInds)-1
            if startInds(jj+1)<endIndsTemp(nextEndInd);
                endInds=[endInds startInds(jj+1)-1];
            else
                nextEndInd=nextEndInd+1;
            end
        end
        
        if length(endInds)<length(startInds)
            endInds=[endInds size(velFolded,1)];
        end
        
        % Go through indices
        for jj=1:length(startInds)
            velVecOut(startInds(jj):endInds(jj))=velVec(startInds(jj):endInds(jj))-2*foldThresh;
        end
        if goodCols(ii)>1
            velCheck=median(velVecOut-velFolded(:,goodCols(ii)-1),'omitnan');
            
            if abs(velCheck)>10
                for jj=1:length(startInds)
                    velVecOut(startInds(jj):endInds(jj))=velVec(startInds(jj):endInds(jj))+2*foldThresh;
                end
                velVecOut=velVecOut-2*foldThresh;
            end
        end
        velUnfold(:,goodCols(ii))=velVecOut;
    end
end
end

