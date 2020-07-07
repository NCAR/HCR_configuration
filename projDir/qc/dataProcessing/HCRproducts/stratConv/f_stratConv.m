% Find stratiform and convective
function [stratConv liquidAlt]= f_meltLayer_altOnly(data,puzzle,meltLayer)

disp('Finding stratiform and convective echo ...');

% Initialize output
stratConvC=nan(size(data.DBZ));
stratConv=nan(size(data.DBZ));

% Find melting layer altitude
meltInds=find(meltLayer>0);
meltAltIn=data.asl(meltInds);

[meltInR meltInC]=ind2sub(size(data.DBZ),meltInds);

meltAlt=nan(size(data.time))';
meltAlt(:,:)=100000;
meltAlt(meltInC)=meltAltIn;

meltR=nan(size(data.time))';
meltR(meltInC)=meltInR;

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
    
    % Stratiform because melting layer detected, but not if cloud is
    % below melting layer
    meltYes=any(meltLayer==1,1);
    meltCloud=nan(size(data.time));
    meltCloud(cloudCu)=meltYes(cloudCu);
    stratConv1D(meltCloud==1 & isnan(stratConv1D))=0;
        
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
    largeFallMask(velMask<=-4)=1;
    
    largeFall=sum(largeFallMask,1);
    stratConv1D(largeFall>3)=1;
    
    % Sea surface cals and antenna transition are unknown
    meanAntFlag=movmean(data.ANTFLAG,500);
    stratConv1D(meanAntFlag>2.5 & (meltAlt'<=maxAltCloud))=2;
    
    % Make it 2D
    backMask=repmat(stratConv1D,size(aslMask,1),1);
    
    % Stratiform because majority of data above melting layer
    maskMask=zeros(size(aslMask));
    maskMask(~isnan(aslMask))=1;
    maskMask=imfill(maskMask,'holes');
        
    for jj=1:size(maskMask,2)
        mRay=maskMask(:,jj);
        altRay=aslMask(:,jj);
        if data.elevation(jj)>0;
            mRay=flipud(mRay);
            altRay=flipud(altRay);
        end
        if ~isnan(meltR(jj))
            mRay(1:meltR(jj))=0;
            firstNonNanM=min(find(mRay==1));
            mRay(1:firstNonNanM-1)=1;
            firstNanM=min(find(mRay==0));
            if ~isempty(firstNanM) & firstNanM~=1 & altRay(firstNanM-1)>minAltCloud(jj)
                newMinAlt=altRay(firstNanM-1);
            else
                newMinAlt=minAltCloud(jj);
            end
        else
            newMinAlt=minAltCloud(jj);
        end
        aboveMelt=maxAltCloud(jj)-meltAlt(jj);
        totCloud=maxAltCloud(jj)-newMinAlt;
        altRayReal=aslMask(:,jj);
        if (aboveMelt./totCloud)>0.8
            backMask(find(altRayReal>=newMinAlt),jj)=0;
            backMask(find(altRayReal<newMinAlt),jj)=4;
        end
    end
            
    % Data above melt (strat) and below melt but not at melt
    noDataMeltIn=isnan(aslMask(meltInds));
    noDataMelt=nan(size(data.time))';
    noDataMelt(:,:)=100000;
    noDataMelt(meltInC)=noDataMeltIn;
    noMeltCols=find(noDataMelt' & meltAlt'>minAltCloud & meltAlt'<maxAltCloud);
    
    for jj=1:length(noMeltCols)
        aslRay=aslMask(:,noMeltCols(jj));
        outRay=nan(size(aslRay));
        outRay(aslRay>=meltAlt(noMeltCols(jj)))=0;
        outRay(aslRay<meltAlt(noMeltCols(jj)))=4;
        backMask(:,noMeltCols(jj))=outRay;
    end
    
    stratConvC(~isnan(aslMask))=backMask(~isnan(aslMask));
end

cutSmall=5000;
cutWidth=20;

% Take care of small convective areas
convMask=zeros(size(stratConvC));
convMask(stratConvC==1)=1;

CC = bwconncomp(convMask);

for jj=1:CC.NumObjects
    area=CC.PixelIdxList{jj};
    [aR aC]=ind2sub(size(convMask),area);
    areaWidth=max(aC)-min(aC);
    if length(area)<cutSmall | areaWidth<cutWidth
        rows=unique(aR);
        minColN=max([min(aC)-1,1]);
        maxColN=min([max(aC)+1,size(stratConvC,2)]);
        
        neighbors=[stratConvC(rows,minColN);stratConvC(rows,maxColN)];
        if any(neighbors==0)
            stratConvC(area)=0;
        elseif any(neighbors==2)
            stratConvC(area)=2;
        elseif any(neighbors==4)
            stratConvC(area)=4;
        end
    end
end

% Take care of small stratiform areas
stratMask=zeros(size(stratConvC));
stratMask(stratConvC==0)=1;

CC = bwconncomp(stratMask);

for jj=1:CC.NumObjects
    area=CC.PixelIdxList{jj};
    [aR aC]=ind2sub(size(stratMask),area);
    areaWidth=max(aC)-min(aC);
    if length(area)<cutSmall | areaWidth<cutWidth
        rows=unique(aR);
        minColN=max([min(aC)-1,1]);
        maxColN=min([max(aC)+1,size(stratConvC,2)]);
        
        neighbors=[stratConvC(rows,minColN);stratConvC(rows,maxColN)];
        if any(neighbors==1)
            stratConvC(area)=1;
        elseif any(neighbors==2)
            stratConvC(area)=2;
        elseif any(neighbors==4)
            stratConvC(area)=4;
        end
    end
end

% Take care of small unknown areas
ukMask=zeros(size(stratConvC));
ukMask(stratConvC==2)=1;

CC = bwconncomp(ukMask);

for jj=1:CC.NumObjects
    area=CC.PixelIdxList{jj};
    [aR aC]=ind2sub(size(ukMask),area);
    areaWidth=max(aC)-min(aC);
    if length(area)<cutSmall | areaWidth<cutWidth
        rows=unique(aR);
        minColN=max([min(aC)-1,1]);
        maxColN=min([max(aC)+1,size(stratConvC,2)]);
        
        neighbors=[stratConvC(rows,minColN);stratConvC(rows,maxColN)];
        if any(neighbors==0)
            stratConvC(area)=0;
        elseif any(neighbors==1)
            stratConvC(area)=1;
        elseif any(neighbors==4)
            stratConvC(area)=4;
        end
    end
end

% Take care of number 4
belowMask=zeros(size(stratConvC));
belowMask(stratConvC==4)=1;

CC = bwconncomp(belowMask);

for jj=1:CC.NumObjects
    area=CC.PixelIdxList{jj};
    [aR aC]=ind2sub(size(belowMask),area);
    rows=unique(aR);
    minColN=max([min(aC)-1,1]);
    maxColN=min([max(aC)+1,size(stratConvC,2)]);
    
    neighbors=[stratConvC(rows,minColN);stratConvC(rows,maxColN)];
    stratConvC(area)=ceil(median(neighbors,'omitnan'));
end

stratConvBig=stratConvC;

% Make it vertically consistent
for jj=1:size(stratConvBig,2)
    scRay=stratConvBig(:,jj);
    if data.elevation(jj)<0;
        scRay=flipud(scRay);
    end
    nanMask=zeros(size(scRay));
    nanMask(~isnan(scRay))=1;
    firstNonNan=min(find(~isnan(scRay)));
    firstVal=scRay(firstNonNan);
    nanMask(1:firstNonNan-1)=1;
    firstNan=min(find(nanMask==0));
    scRay(firstNonNan:firstNan)=firstVal;
    if data.elevation(jj)<0;
        scRay=flipud(scRay);
    end
    stratConv(:,jj)=scRay;
end

% Find maximum reflectivity altitude
typeMeltIn=stratConv(meltInds);
typeMelt=nan(size(data.time))';
typeMelt(meltInC)=typeMeltIn;

onlyConv=stratConv;
onlyConv(stratConv~=1)=nan;
onlyConv(onlyConv==1)=data.dbzMasked(onlyConv==1);

[dbzMax dbzMaxInd]=max(onlyConv,[],1);
maxInLin=sub2ind(size(data.DBZ),dbzMaxInd,1:length(data.time));
liquidAlt=data.asl(maxInLin)';

liquidAlt(liquidAlt<=meltAlt+100)=nan; % Remove data that is too close to melting layer
liquidAlt(dbzMaxInd==1)=nan;

% Add data for strat and no data at melting layer
liquidAlt(isnan(liquidAlt) & (typeMelt==0 | isnan(typeMelt)))=meltAlt(isnan(liquidAlt) & (typeMelt==0 | isnan(typeMelt)));

% Remove data when the plane is in the convective cloud
inCloud=stratConv(18,:);
inCloud(data.elevation>0)=nan;
inCloud(inCloud~=1)=nan;

% If plane is really high the data is good
inCloud(data.asl(18,:)>13000)=nan;
liquidAlt(inCloud==1)=nan;

end