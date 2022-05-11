function [zenithCorrection,medCloudVelNadir,medCloudVelZenith,smoothfactor,velZenithCorrSmooth]=velCorrZenith(data,velNadirCorr)

%% Get cloud Top vel
% Mask non-cloud data
velMasked=velNadirCorr;
velMasked(data.FLAG~=1)=nan;
velMasked(:,data.ANTFLAG>2)=nan;

% Reverse sign in zenith pointing
velMasked(:,data.elevation>0)=-velMasked(:,data.elevation>0);

% Find cloud top vel
cloudTopVel=nan(size(velMasked,2),1);
for mm=1:size(velMasked,2)
    if mm==151828
        stophere=1;
    end

    velCol=velMasked(:,mm);
    if all(isnan(velCol))
        continue
    end

    % Find individual clouds
    cloudMask=~isnan(velCol);
    % Join clouds that are close
    cloudMask=movmedian(cloudMask,9,'omitnan');
    % Remove small clouds
    cloudMask=bwareaopen(cloudMask,10);
    if sum(cloudMask)==0
        continue
    end

    % Loop through individual clouds
    cloudTops=[];
    clouds=bwconncomp(cloudMask);
    for nn=1:clouds.NumObjects
        thisCloud=nan(size(cloudMask));
        thisCloud(clouds.PixelIdxList{nn})=velCol(clouds.PixelIdxList{nn});
        
        % Check for many strong updraft pixels
        extVel=sum(thisCloud<-1);
        if extVel/sum(~isnan(thisCloud))>0.1
            continue
        end

        % Collect cloud top data
        if data.elevation(mm)<=0
            % Check if in cloud
            %if ~isnan(thisCloud(22))
            if ~isempty(intersect(clouds.PixelIdxList{nn},22))
                continue
            end
            firstInd=min(find(~isnan(thisCloud)));
            cloudTops=cat(1,cloudTops,thisCloud(firstInd:firstInd+2));
        else
            lastInd=max(find(~isnan(thisCloud)));
            % Check if near surface
            cloudAsl=data.asl(lastInd-2:lastInd,mm);
            if max(cloudAsl,'omitnan')<1000
                continue
            end
            cloudTops=cat(1,cloudTops,thisCloud(lastInd-2:lastInd));
        end
    end
    if isempty(cloudTops)
        continue
    elseif all(~isnan(cloudTops))
        cloudTopVel(mm)=median(cloudTops,1,'omitnan');
    end
end

%% Calc things
smoothfactor=10000;

cloudTopVelInds=find(~isnan(cloudTopVel));
cloudTopNadirInds=find(data.elevation<=0);
cloudTopZenithInds=find(data.elevation>0);

cloudTopVelIndsNadir=intersect(cloudTopVelInds,cloudTopNadirInds);
cloudTopVelShrunkNadir=cloudTopVel(cloudTopVelIndsNadir);

medCloudVelShrunkNadir=movmedian(cloudTopVelShrunkNadir,smoothfactor,'omitnan');
medCloudVelNadir=nan(size(cloudTopVel));
medCloudVelNadir(cloudTopVelIndsNadir)=medCloudVelShrunkNadir;
medCloudVelNadirInt=fillmissing(medCloudVelNadir,'linear','EndValues','nearest');

cloudTopVelIndsZenith=intersect(cloudTopVelInds,cloudTopZenithInds);
cloudTopVelOnlyZenith=nan(size(cloudTopVel));
cloudTopVelOnlyZenith(cloudTopVelIndsZenith)=cloudTopVel(cloudTopVelIndsZenith);

% Remove isolated
ctvelz=movmean(cloudTopVelOnlyZenith,100,'omitnan');
ctMask=~isnan(ctvelz);
ctMask=bwareaopen(ctMask,180);
cloudTopVelOnlyZenith(ctMask==0)=nan;

medCloudVelZenith=movmedian(cloudTopVelOnlyZenith,smoothfactor*2,'omitnan');
medCloudVelZenith(isnan(cloudTopVelOnlyZenith))=nan;

zenithCorrection=medCloudVelNadirInt-medCloudVelZenith;
zenithCorrection=fillmissing(zenithCorrection,'linear','EndValues','nearest');
zenithCorrection(data.elevation<=0)=nan;
velZenithCorrSmooth=movmedian(cloudTopVelOnlyZenith,smoothfactor,'omitnan')+zenithCorrection;
end