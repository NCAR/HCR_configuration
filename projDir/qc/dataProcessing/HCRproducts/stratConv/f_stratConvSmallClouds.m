function [stratConv,stratConv1D]=f_stratConvSmallClouds(stratConv,stratConv1D,cloudPuzzle,dbzSmallClouds,asl,icingLevel)
% Stratiform/convective partitioning for small clouds that are marked as 0
% in cloud puzzle
zeroClouds=zeros(size(stratConv));
zeroClouds(cloudPuzzle==0)=1;
asl(zeroClouds==0)=nan;

indivClouds=bwconncomp(zeroClouds);

% Cloud pixel threshold for really small clouds
reallySmall=1000;

% Loop through clouds
for ii=1:indivClouds.NumObjects
    cloudPix=indivClouds.PixelIdxList{ii};
    % If cloud is really small classify as convective (stratiform) below
    % (above) the icing level
    [cloudR,cloudC]=ind2sub(size(stratConv),cloudPix);
    if length(cloudPix)<reallySmall
        % Count number of pixels above and below icing level
        uniqueCols=unique(cloudC);
        aboveI=0;
        belowI=0;
        for jj=1:length(uniqueCols)
            altCol=asl(:,uniqueCols(jj));
            altCol(isnan(altCol))=[];
            aboveI=aboveI+sum(altCol>icingLevel(uniqueCols(jj)));
            belowI=belowI+sum(altCol<icingLevel(uniqueCols(jj)));
        end
        if aboveI>=belowI
            stratConv(cloudPix)=2;
        else
            stratConv(cloudPix)=1;
        end
    else
        % If height>width -> convective
        heightI=abs(max(cloudR)-min(cloudR));
        widthI=max(cloudC)-min(cloudC);
        if heightI>widthI
            stratConv(cloudPix)=1;
        else
            stratConv(cloudPix)=2;
        end
    end
    % If more than x% of dbz data is >x dBZ -> convective
    reflCut=5;
    reflPerc=0.03;
    if length(find(dbzSmallClouds(cloudPix)>reflCut))/length(cloudPix)>reflPerc
        stratConv(cloudPix)=1;
    end
end

% Add new 1D values where cloud classifications were added
% Convert 2D stratiform convective field to 1D
stratConv1Dnew=nan(1,size(stratConv,2));
stratConv1Dnew(:)=2;

stratConvBelow=stratConv;
stratCount=zeros(size(stratConv1Dnew));
for ii=1:length(stratConv1Dnew)
    aslCol=asl(:,ii);
    stratConvBelow(aslCol>icingLevel(ii),ii)=nan;
    
    stratCount(ii)=sum(stratConv(:,ii)==2);
end

stratConv1Dnew(any(stratConvBelow==1,1))=1;
stratConv1Dnew(stratConv1Dnew==2 & stratCount==0)=1;

zeroClouds=bwareaopen(zeroClouds,reallySmall);
stratConv1Dnew(sum(zeroClouds,1)==0)=nan;

stratConv1D(isnan(stratConv1D))=stratConv1Dnew(isnan(stratConv1D));
end