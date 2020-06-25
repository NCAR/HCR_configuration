% Find stratiform and convective
function stratConv= f_meltLayer_altOnly(data,puzzle,meltLayer,meltArea)

disp('Finding stratiform and convective echo ...');

% Initialize output
stratConv=nan(size(data.DBZ));

% Find melting layer altitude
meltInds=find(meltLayer>0);
meltAlt=data.asl(meltInds);

[meltInR meltInC]=ind2sub(size(data.DBZ),meltInds);

% Loop through clouds
numClouds=max(max(puzzle));
countPieces=1;
newPuzzle=nan(size(puzzle));
keepCols=[];

% Find cloud layers and split puzzle even further
for ii=1:numClouds
    cloudInd=find(puzzle==ii);
    
    cloudMask=zeros(size(puzzle));
    cloudMask(cloudInd)=1;
    
    % Fill holes
    cloudMaskFilled= imfill(cloudMask,'holes');
    
    for jj=1:size(cloudMaskFilled,2)
        rayInds=find(cloudMaskFilled(:,jj)==1);
        min1=min(rayInds);
        max1=max(rayInds);
        
        cloudMaskFilled(1:min1,jj)=1;
        cloudMaskFilled(max1:end,jj)=1;
    end
    
    % Reverse image and remove small objects
    cloudMaskRev=~cloudMaskFilled;
    cloudMaskBig=bwareaopen(cloudMaskRev,5000);
    
    % Loop through contiguous areas
    CC = bwconncomp(cloudMaskBig);
        
    for jj=1:CC.NumObjects
        area=CC.PixelIdxList{jj};
        [areaRows areaCols]=ind2sub(size(cloudMaskBig),area);
                
        % Width of piece
        layerWidth=max(areaCols)-min(areaCols);
        if layerWidth>500
            cloudMask(:,min(areaCols))=0;
            cloudMask(:,max(areaCols))=0;
            %keepCols=[keepCols,min(areaCols),max(areaCols)];
        end
    end
    
    % Put new pieces in new puzzle
    CC2 = bwconncomp(cloudMask);
        
    for jj=1:CC2.NumObjects
        area2=CC2.PixelIdxList{jj};
        newPuzzle(area2)=countPieces;
        
        % Fill in the empty column
        [areaRows2 areaCols2]=ind2sub(size(cloudMask),area2);
        
        maxCol=max(areaCols2);
        if maxCol~=size(cloudMask,2);
            colData=newPuzzle(:,maxCol);
            rowInds=find(colData==countPieces);
            newPuzzle(rowInds,maxCol+1)=countPieces;
        end
        
        countPieces=countPieces+1;
    end
end

for ii=1:countPieces-1
    stratConv1D=nan(size(data.time));
    
    cloudInd=find(newPuzzle==ii);
    [cloudR cloudC]=ind2sub(size(data.dbzMasked),cloudInd);
    cloudCu=unique(cloudC);
    cloudYes=zeros(size(data.time));
    cloudYes(cloudCu)=1;
    
    aslMask=nan(size(data.dbzMasked));
    aslMask(cloudInd)=data.asl(cloudInd);
    
    maxAltCloud=max(aslMask,[],1);
    minAltCloud=min(aslMask,[],1);
    
    % Shallow convection
    stratConv1D(meltAlt'>maxAltCloud)=1;
    
    % Stratiform because melting layer detected
    meltYes=any(meltLayer==1,1);
    meltCloud=nan(size(data.time));
    meltCloud(cloudCu)=meltYes(cloudCu);
    stratConv1D(meltCloud==1)=0;
        
    % Convective, no melting layer found
    stratConv1D(cloudYes & isnan(stratConv1D))=1;
    
    % But if VEL above melting layer > 2.5 m/s -> convective
    velMask=nan(size(data.dbzMasked));
    velMask(cloudInd)=data.VEL_CORR(cloudInd);
    
    for jj=1:size(velMask,2);
        aslRay=aslMask(:,jj);
        velMask(aslRay<=meltAlt(jj)+200,jj)=nan;
    end
    
    velMask(:,data.elevation<0)=-velMask(:,data.elevation<0);   
    largeFallMask=zeros(size(velMask));
    largeFallMask(velMask<=-2)=1;
    
    largeFall=sum(largeFallMask,1);
    stratConv1D(largeFall>3)=1;
    
    % Stratiform because above melting layer
    stratConv1D(meltAlt'<minAltCloud)=0;   
    
    % Make it 2D
    backMask=repmat(stratConv1D,size(aslMask,1),1);
    
    % Data above melt and below melt but not at melt
    noDataMelt=isnan(aslMask(meltInds));
    noMeltCols=find(noDataMelt' & meltAlt'>minAltCloud & meltAlt'<maxAltCloud);
        
    for jj=1:length(noMeltCols)
        aslRay=aslMask(:,noMeltCols(jj));
        outRay=nan(size(aslRay));
        outRay(aslRay>=meltAlt(noMeltCols(jj)))=0;
        outRay(aslRay<meltAlt(noMeltCols(jj)))=2;
        backMask(:,noMeltCols(jj))=outRay;
    end
    
    % Get rid of outliers
    backMask(~isnan(aslMask))=backMask(~isnan(aslMask));
    stratConvMed=floor(movmedian(backMask,49,2,'omitnan'));
    stratConvMed=floor(movmedian(stratConvMed,9,1,'omitnan'));
    
    stratConv(~isnan(aslMask))=stratConvMed(~isnan(aslMask));
end

% Take care of number 2
belowMask=zeros(size(stratConv));
belowMask(stratConv==2)=1;

CC = bwconncomp(belowMask);

for jj=1:CC.NumObjects
    area=CC.PixelIdxList{jj};
    [aR aC]=ind2sub(size(belowMask),area);
    rows=unique(aR);
    minColN=max([min(aC)-1,1]);
    maxColN=min([max(aC)+1,size(stratConv,2)]);
    
    neighbors=[stratConv(rows,minColN);stratConv(rows,maxColN)];
    stratConv(area)=ceil(median(neighbors,'omitnan'));
end

% % Find maximum reflectivity altitude
% [dbzMax dbzMaxInd]=max(dbzSmooth,[],1);
% maxInLin=sub2ind(size(data.DBZ),dbzMaxInd,1:length(data.time));
% maxAlt=data.asl(maxInLin);
% 
% % Find melting layer altitude
% meltInds=find(meltLayer>0);
% meltAlt=data.asl(meltInds);
% 
% % Distance between max refl and melt alt
% reflDist=maxAlt-meltAlt';
% stratConv(reflDist>500)=1;
% stratConv(reflDist<=500)=0;
% 
% % Is there data at melting layer?
% meltData=data.dbzMasked(meltInds);
% 
% stratConv(isnan(meltData') & reflDist>0)=0;
% stratConv(isnan(meltData') & reflDist<0)=2;

end