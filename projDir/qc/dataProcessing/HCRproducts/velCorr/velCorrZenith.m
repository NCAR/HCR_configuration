function [zenithCorrection,medCloudVelNadir,medCloudVelZenith,smoothfactor,velZenithCorrSmooth]=velCorrZenith(data,velNadirCorr)

%% Get cloud Top vel
% Mask non-cloud data
velMasked=velNadirCorr;
velMasked(data.FLAG~=1)=nan;
velMasked(:,data.ANTFLAG>2)=nan;

% Reverse sign in zenith pointing
velMasked(:,data.elevation>0)=-velMasked(:,data.elevation>0);

% Remove small clouds
maskVel=~isnan(velMasked);
maskVel=bwareaopen(maskVel,1000);

% Join clouds that are close
maskVel=imdilate(maskVel,strel('disk',9));
labs=bwlabel(maskVel);
labs(labs==0)=nan;
labs(isnan(velMasked))=nan;

% Find cloud top vel
cloudTopVel=nan(size(velMasked));
cloudTopAlts=nan(size(velMasked));
for mm=1:size(velMasked,2)
    if mm==116155
        stophere=1;
    end

    labCol=labs(:,mm);
    if all(isnan(labCol))
        continue
    end
    velCol=velMasked(:,mm);
    altCol=data.asl(:,mm);
    
    % Find individual clouds
    cloudMask=~isnan(labCol);
    % Remove small clouds
    cloudMask=bwareaopen(cloudMask,10);
    if sum(cloudMask)==0
        continue
    end
    labCol(cloudMask==0)=nan;
    velCol(cloudMask==0)=nan;

    % Loop through individual clouds
%     cloudTops=[];
%     cloudAlts=[];
    clouds=unique(labCol);
    clouds(isnan(clouds))=[];
    %clouds=bwconncomp(cloudMask);
    for nn=1:length(clouds)
        thisCloud=nan(size(cloudMask));
        thisCloud(labCol==clouds(nn))=velCol(labCol==clouds(nn));
        
        % Check for many strong updraft pixels
        extVel=sum(thisCloud<-1);
        if extVel/sum(~isnan(thisCloud))>0.1
            continue
        end

        % Collect cloud top data
        if data.elevation(mm)<=0
            % Check if in cloud
            if ~isempty(intersect(find(labCol==clouds(nn)),22))
                continue
            end
            firstIndFirst=min(find(~isnan(thisCloud)));
            thisCloud(abs(thisCloud)>2)=nan;
            newMask=~isnan(thisCloud);
            newMask=bwareaopen(newMask,10);
            thisCloud(newMask==0)=nan;
            firstInd=min(find(~isnan(thisCloud)));
            if abs(firstIndFirst-firstInd)>15
                continue
            end
            cloudTopVel(firstInd,mm)=median(thisCloud(firstInd:firstInd+2),1,'omitnan');
            cloudTopAlts(firstInd,mm)=data.asl(firstInd,mm);
%             cloudTops=cat(1,cloudTops,thisCloud(firstInd:firstInd+2));
%             cloudAlts=cat(1,cloudAlts,altCol(firstInd:firstInd+2));
        else
            lastIndFirst=max(find(~isnan(thisCloud)));
            thisCloud(abs(thisCloud)>2)=nan;
            newMask=~isnan(thisCloud);
            newMask=bwareaopen(newMask,10);
            thisCloud(newMask==0)=nan;
            lastInd=max(find(~isnan(thisCloud)));
            if abs(lastIndFirst-lastInd)>15
                continue
            end
            % Check if near surface
            cloudAsl=data.asl(lastInd-2:lastInd,mm);
            if max(cloudAsl,'omitnan')<1000
                continue
            end
            cloudTopVel(lastInd,mm)=median(thisCloud(lastInd-2:lastInd),1,'omitnan');
            cloudTopAlts(lastInd,mm)=data.asl(lastInd,mm);
%             cloudTops=cat(1,cloudTops,thisCloud(lastInd-2:lastInd));
%             cloudAlts=cat(1,cloudAlts,altCol(lastInd-2:lastInd));
        end
    end
%     if isempty(cloudTops)
%         continue
%     elseif all(~isnan(cloudTops))
%         altRound=round(cloudAlts./1000).*1000;
%         cloudTopAlts(mm)=mode(altRound);
%         cloudTops(abs(altRound-cloudTopAlts(mm))>3000)=[];
%         cloudTopVel(mm)=median(cloudTops,1,'omitnan');
%     end
end

%% Loop through clouds

ulabs=unique(labs);
ulabs(isnan(ulabs))=[];

allAsls=data.asl(~isnan(cloudTopAlts));
velLayers=nan(ceil(max(allAsls)/1000),size(velMasked,2));

for ii=1:length(ulabs)
    labInds=find(labs==ulabs(ii));
    velTopLabInds=cloudTopVel(labInds);
    labInds(isnan(velTopLabInds))=[];
    velTopLabInds(isnan(velTopLabInds))=[];
    if isempty(velTopLabInds)
        continue
    end
    medAltCl=median(data.asl(labInds),'omitnan');
    medAltRound=round(medAltCl/1000);

    [labR,labC]=ind2sub(size(velMasked),labInds);

    velLayers(medAltRound,labC)=velTopLabInds;
end

nadirInds=find(data.elevation<0);
zenithInds=find(data.elevation>=0);
% 
% nonNanInds=find(any(~isnan(velLayers),1));
% nadirCloudInds=intersect(nadirInds,nonNanInds);
% zenithCloudInds=intersect(zenithInds,nonNanInds);

velDown=velLayers;
velDown(:,zenithInds)=nan;
velUp=velLayers;
velUp(:,nadirInds)=nan;

smoothfactor=10001;
se=strel('rectangle',[1,smoothfactor]);

velDownMask=~isnan(velDown);
velDownMask=imdilate(velDownMask,se);
connDown=bwconncomp(velDownMask);
velDownGood=nan(size(velDown));

for ii=1:connDown.NumObjects
    thisReg=nan(size(velLayers));
    thisReg(connDown.PixelIdxList{ii})=velDown(connDown.PixelIdxList{ii});
    thisRegVel=median(thisReg,1,'omitnan');
    thisVelSmooth=movmedian(thisRegVel,smoothfactor,'omitnan');
    thisVelSmooth(isnan(thisRegVel))=nan;
    thisVelInt=fillmissing(thisVelSmooth,'linear','EndValues','nearest');
    [rthis,cthis]=ind2sub(size(velLayers),connDown.PixelIdxList{ii});
    uc=unique(cthis);
    for jj=1:length(uc)
        velDownGood(find(velDownMask(:,uc(jj))>0),uc(jj))=thisVelInt(uc(jj));
    end
end

velUpMask=~isnan(velUp);
velUpMask=imdilate(velUpMask,se);
connUp=bwconncomp(velUpMask);
%velUpGood=nan(size(velUp));

for ii=1:connUp.NumObjects
    thisReg=nan(size(velLayers));
    thisReg(connUp.PixelIdxList{ii})=velUp(connUp.PixelIdxList{ii});
    thisRegVel=median(thisReg,1,'omitnan');
    thisVelSmooth=movmedian(thisRegVel,smoothfactor,'omitnan');
    thisVelSmooth(isnan(thisRegVel))=nan;
    thisVelInt=fillmissing(thisVelSmooth,'linear','EndValues','nearest');
    [rthis,cthis]=ind2sub(size(velLayers),connUp.PixelIdxList{ii});

    % Find closest nadir region
    minRow=min(rthis);
    maxRow=max(rthis);

%     uc=unique(cthis);
%     for jj=1:length(uc)
%         velUpGood(find(velUpMask(:,uc(jj))>0),uc(jj))=thisVelInt(uc(jj));
%     end
end


%% Make sure we are at the right altitude
smoothfactor=10001;
cloudTopRemove=cloudTopAlts;
cloudTopRemove(isnan(cloudTopVel))=nan;
cloudTopRemove(data.elevation>0)=nan;
movAlt=movmedian(cloudTopRemove,smoothfactor,'omitnan');
% movAlt(isnan(cloudTopVel))=nan;
% movAlt(data.elevation>0)=nan;
movAlt=fillmissing(movAlt,'linear','EndValues','nearest');

%% Calc things

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
ctvelz=movmean(cloudTopVelOnlyZenith,101,'omitnan');
ctMask=~isnan(ctvelz);
ctMask=bwareaopen(ctMask,180);
cloudTopVelOnlyZenith(ctMask==0)=nan;

% Remove those that were collected at the wrong altitude
cloudTopVelOnlyZenith(abs(cloudTopAlts-movAlt)>3000)=nan;

medCloudVelZenith=movmedian(cloudTopVelOnlyZenith,smoothfactor*2+1,'omitnan');
medCloudVelZenith(isnan(cloudTopVelOnlyZenith))=nan;

zenithCorrection=medCloudVelNadirInt-medCloudVelZenith;
zenithCorrection=fillmissing(zenithCorrection,'linear','EndValues','nearest');
zenithCorrection(data.elevation<=0)=nan;
velZenithCorrSmooth=movmedian(cloudTopVelOnlyZenith,smoothfactor,'omitnan')+zenithCorrection;
end