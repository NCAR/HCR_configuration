function velDeAlias=dealiasArea(velFolded,elev,nyq)
% Unfold velocities

% Make all updrafts positive and downdrafts negative
velFolded(:,elev<0)=-velFolded(:,elev<0);
velDeAlias=nan(size(velFolded));

% Check if nyquist velocity changes
if max(nyq)-min(nyq)~=0
    error('Nyquist velocity is not constant.')
end

% Go through areas
velMask=~isnan(velFolded);
velAreas=bwconncomp(velMask);

for ii=1:velAreas.NumObjects
    
    if rem(ii,50)==0
        disp(['De-aliasing area ',num2str(ii),' of ',num2str(velAreas.NumObjects),'.']);
    end
    
    % Cut out small area
    [rowInds colInds]=ind2sub(size(velMask),velAreas.PixelIdxList{ii});
    oneVel=nan(size(velFolded));
    oneVel(velAreas.PixelIdxList{ii})=velFolded(velAreas.PixelIdxList{ii});
    oneVel=oneVel(min(rowInds):max(rowInds),min(colInds):max(colInds));
        
    % Dealias in the positive direction
    [velDeAliasOne doNeg]=dealiasAreaPos(oneVel,nyq(1));
    
    if doNeg
        velDeAliasOne=dealiasAreaNeg(velDeAliasOne,nyq(1));
        [velDeAliasOne ~]=dealiasAreaPos(velDeAliasOne,nyq(1));
    end
    
    velDeAlias(velAreas.PixelIdxList{ii})=velDeAliasOne(~isnan(velDeAliasOne));
end

% Transform down pointing back to original sign
velDeAlias(:,elev<0)=-velDeAlias(:,elev<0);

end