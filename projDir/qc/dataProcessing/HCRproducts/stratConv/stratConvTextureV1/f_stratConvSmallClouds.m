function [stratConv,stratConv1D]=f_stratConvSmallClouds(stratConv,stratConv1D,cloudPuzzle,dbzSmallClouds,meltingLayer)
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
    
    % Count number of pixels above and below icing level
    meltCloud=meltingLayer(cloudPix);
    aboveI=sum(meltCloud==20);
    belowI=sum(meltCloud==10);
    
    % If cloud is really small classify as convective (stratiform) below
    % (above) the icing level
    if length(cloudPix)<reallySmall
        if aboveI>=belowI
            stratConv(cloudPix)=22;
        else
            stratConv(cloudPix)=13;
        end
    else
        % If height>width -> convective
        [cloudR,cloudC]=ind2sub(size(stratConv),cloudPix);
        heightI=abs(max(cloudR)-min(cloudR));
        widthI=max(cloudC)-min(cloudC);
        if heightI>widthI/2
            if aboveI>=belowI
                stratConv(cloudPix)=14;
            else
                stratConv(cloudPix)=13;
            end
        else
            if aboveI>=belowI
                stratConv(cloudPix)=22;
            else
                stratConv(cloudPix)=23;
            end
        end
    end
    % If more than x% of dbz data is >x dBZ -> convective
    reflCut=5;
    reflPerc=0.03;
    if length(find(dbzSmallClouds(cloudPix)>reflCut))/length(cloudPix)>reflPerc
        if aboveI>=belowI
            stratConv(cloudPix)=14;
        else
            stratConv(cloudPix)=13;
        end
    end
end

% Add new 1D values where cloud classifications were added
% Convert 2D stratiform convective field to 1D
stratConv1Dnew=nan(1,size(stratConv,2));
stratConv1Dnew(:)=2;

zeroClouds2=bwareaopen(zeroClouds,reallySmall);
stratConv2=stratConv;
stratConv2(zeroClouds==1 & zeroClouds2==0)=nan;

stratConvBelow=stratConv2;
stratConvBelow(meltingLayer~=10)=nan;
stratCount=sum(stratConv2>21,1);

stratConv1Dnew(any(stratConvBelow==13,1))=1;
stratConv1Dnew(stratConv1Dnew==2 & stratCount==0)=1;

stratConv1Dnew(sum(zeroClouds2,1)==0)=nan;

stratConv1D(isnan(stratConv1D))=stratConv1Dnew(isnan(stratConv1D));
end