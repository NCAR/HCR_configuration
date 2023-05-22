function [zenithCorrection,smoothfactor]=velCorrZenith(data,velNadirCorr)

%% Clean up input data
% Mask non-cloud data
velMasked=velNadirCorr;
velMasked(data.FLAG~=1)=nan;
velMasked(:,data.ANTFLAG>2)=nan;

% Reverse sign in zenith pointing
velMasked(:,data.elevation>0)=-velMasked(:,data.elevation>0);

% Remove small clouds
maskVel=~isnan(velMasked);
maskVel=bwareaopen(maskVel,5000);
velMasked(maskVel==0)=nan;

% Join clouds that are close
maskVel=imdilate(maskVel,strel('disk',9));

% Create label mat with individual clouds
labs=bwlabel(maskVel);
labs(labs==0)=nan;
labs(isnan(velMasked))=nan;

%% Find cloud top vel
% 2D mat with cloud top velocities
cloudTopVel=nan(size(velMasked));
% 2D mat with altitudes of cloud tops
cloudTopAlts=nan(size(velMasked));
for mm=1:size(velMasked,2)
    if mm==137342
        stophere=1;
    end

    % Create column vectors
    labCol=labs(:,mm);
    if all(isnan(labCol))
        continue
    end
    velCol=velMasked(:,mm);
        
    % Find individual cloud pieces in col vectors
    cloudMask=~isnan(labCol);
    % Remove small cloud pieces
    cloudMask=bwareaopen(cloudMask,10);
    if sum(cloudMask)==0
        continue
    end
    labCol(cloudMask==0)=nan;
    velCol(cloudMask==0)=nan;

    clouds=unique(labCol);
    clouds(isnan(clouds))=[];
    
    % Loop through cloud pieces to get top vel and altitude of top
    for nn=1:length(clouds)
        thisCloud=nan(size(cloudMask));
        thisCloud(labCol==clouds(nn))=velCol(labCol==clouds(nn));
        
        % Check for many strong updraft pixels and exclude cloud if found
        extVel=sum(thisCloud<-1);
        if extVel/sum(~isnan(thisCloud))>0.1
            continue
        end

        % Collect cloud top data
        if data.elevation(mm)<=0 % nadir
            % Check if in cloud
            if ~isempty(intersect(find(labCol==clouds(nn)),22))
                continue
            end
            firstIndFirst=min(find(~isnan(thisCloud)));
            % Remove noise pixels (i.e. pixels with high vels at cloud edge)
            thisCloud(abs(thisCloud)>2)=nan;
            newMask=~isnan(thisCloud);
            newMask=bwareaopen(newMask,10);
            thisCloud(newMask==0)=nan;
            % Top index
            firstInd=min(find(~isnan(thisCloud)));
            % Check if noise remove step removed pixels too far away from the edge
            if abs(firstIndFirst-firstInd)>15
                continue
            end
            % Take vels of top three
            cloudTopVel(firstInd,mm)=median(thisCloud(firstInd:firstInd+2),1,'omitnan');
            cloudTopAlts(firstInd,mm)=data.asl(firstInd,mm);
        else % zenith
            lastIndFirst=max(find(~isnan(thisCloud)));
            % Remove noise pixels (i.e. pixels with high vels at cloud edge)
            thisCloud(abs(thisCloud)>2)=nan;
            newMask=~isnan(thisCloud);
            newMask=bwareaopen(newMask,10);
            thisCloud(newMask==0)=nan;
            % Top index
            lastInd=max(find(~isnan(thisCloud)));
            % Check if noise remove step removed pixels too far away from the edge
            if abs(lastIndFirst-lastInd)>15
                continue
            end
            % Check if near surface
            cloudAsl=data.asl(lastInd-2:lastInd,mm);
            if max(cloudAsl,'omitnan')<500
                continue
            end
            % Take vels of top three
            cloudTopVel(lastInd,mm)=median(thisCloud(lastInd-2:lastInd),1,'omitnan');
            cloudTopAlts(lastInd,mm)=data.asl(lastInd,mm);
        end
    end
end

%% Collect cloud top vels and altitudes in altitude bins of 1 km

ulabs=unique(labs);
ulabs(isnan(ulabs))=[];

% Absolute maximum altitude
allAsls=data.asl(~isnan(cloudTopAlts));
% Initiate altitude bin mat
velLayers=nan(ceil(max(allAsls)/1000),size(velMasked,2));

% Loop through clouds and put cloud top vels in altitude bin mat
for ii=1:length(ulabs)
    % Indices of cloud
    labInds=find(labs==ulabs(ii));
    % Remove non cloud top indices
    velTopLabInds=cloudTopVel(labInds);
    labInds(isnan(velTopLabInds))=[];
    velTopLabInds(isnan(velTopLabInds))=[];
    if isempty(velTopLabInds)
        continue
    end
    % Find altitude bin by calculating median altitude of cloud tops
    medAltCl=median(data.asl(labInds),'omitnan');
    medAltRound=max([round(medAltCl/1000),1]);

    % Put cloud top vel in correct altitude bin
    [~,labC]=ind2sub(size(velMasked),labInds);
    velLayers(medAltRound,labC)=velTopLabInds;
end

%% Separate altitude bin vels in nadir and zenith
velDown=velLayers;
velDown(:,data.elevation>=0)=nan;
velUp=velLayers;
velUp(:,data.elevation<0)=nan;

%% Connect altitudes in alt bin mat that are close together (nadir)

% Spread the alt bin mat data out mostly horizontally, but also vertically,
% so close-by regions connect
smoothfactor=round(hcrTimeToPix(30,etime(datevec(data.time(2)),datevec(data.time(1))))+1);
se=strel('rectangle',[1,round(smoothfactor/2)]);

velDownMask=~isnan(velDown);
velDownMask=imdilate(velDownMask,se);
connDown=bwconncomp(velDownMask);

% Find separate levels of cloud tops
velDownLayers=nan(1,size(velDown,2));
altDownLayers=nan(1,size(velDown,2));

% Loop through clouds
for ii=1:connDown.NumObjects
    thisRegVel=nan(size(velLayers));
    thisRegVel(connDown.PixelIdxList{ii})=velDown(connDown.PixelIdxList{ii});
    % Vertical mean over top three pixels
    thisRegVel=median(thisRegVel,1,'omitnan');
    % Smooth
    thisVelSmooth=movmedian(thisRegVel,smoothfactor,'omitnan');
    thisVelSmooth(isnan(thisRegVel))=nan;

    % Find altitude
    [rthis,~]=ind2sub(size(velLayers),connDown.PixelIdxList{ii});
    meanAlt=mean(rthis); % Altitude of this cloud
    if ii==1
        altDownLayers(~isnan(thisVelSmooth))=meanAlt;
        velDownLayers=thisVelSmooth;
    else
        % Compare with existing cloud levels
        diffAlts=abs(allMeanAltDown-meanAlt);
        [minDiffAlt,minInd]=min(diffAlts);
        if minDiffAlt<=3 % If altitude differenc is small, add to existing level
            altDownLayers(minInd,~isnan(thisVelSmooth))=meanAlt;
            velDownLayers(minInd,~isnan(thisVelSmooth))=thisVelSmooth(~isnan(thisVelSmooth));
        else % Otherwise create new level
            altAddRow=nan(1,size(velDown,2));
            altAddRow(~isnan(thisVelSmooth))=meanAlt;
            altDownLayers=cat(1,altDownLayers,altAddRow);
            velDownLayers=cat(1,velDownLayers,thisVelSmooth);
        end
    end
    allMeanAltDown=mean(altDownLayers,2,'omitnan');
end

%% Repeat with zenith

% Spread the alt bin mat data out mostly horizontally, but also vertically,
% so close-by regions connect
velUpMask=~isnan(velUp);
velUpMask=imdilate(velUpMask,se);
connUp=bwconncomp(velUpMask);

if connUp.NumObjects==0
    zenithCorrection=zeros(1,size(velNadirCorr,2));
    disp('No zenith correction found.');
    return
end

% Find separate levels of cloud tops
velUpLayers=nan(1,size(velUp,2));
altUpLayers=nan(1,size(velUp,2));

for ii=1:connUp.NumObjects
    thisRegVel=nan(size(velLayers));
    thisRegVel(connUp.PixelIdxList{ii})=velUp(connUp.PixelIdxList{ii});
    % Vertical mean over top three pixels
    thisRegVel=median(thisRegVel,1,'omitnan');
    % Smooth
    thisVelSmooth=movmedian(thisRegVel,smoothfactor,'omitnan');
    thisVelSmooth(isnan(thisRegVel))=nan;

    % Find altitude
    [rthis,cthis]=ind2sub(size(velLayers),connUp.PixelIdxList{ii});
    meanAlt=mean(rthis);
    if ii==1
        altUpLayers(~isnan(thisVelSmooth))=meanAlt;
        velUpLayers=thisVelSmooth;
    else % Compare with existing cloud levels
        diffAlts=abs(allMeanAltUp-meanAlt);
        [minDiffAlt,minInd]=min(diffAlts);
        if minDiffAlt<=3
            altUpLayers(minInd,~isnan(thisVelSmooth))=meanAlt;
            velUpLayers(minInd,~isnan(thisVelSmooth))=thisVelSmooth(~isnan(thisVelSmooth));
        else % Otherwise create new level
            altAddRow=nan(1,size(velUp,2));
            altAddRow(~isnan(thisVelSmooth))=meanAlt;
            altUpLayers=cat(1,altUpLayers,altAddRow);
            velUpLayers=cat(1,velUpLayers,thisVelSmooth);
        end
    end
    allMeanAltUp=mean(altUpLayers,2,'omitnan');
end

%% Create correction
% Remove small
maskDown=movmean(velDownLayers,5,2,'omitnan');
maskDown=movmean(maskDown,25,2,'includenan');
velDownLayers(isnan(maskDown))=nan;
maskUp=movmean(velUpLayers,5,2,'omitnan');
maskUp=movmean(maskUp,25,2,'includenan');
velUpLayers(isnan(maskUp))=nan;

% Connect cloud top vels horizontally in the different levels
downInt=fillmissing(velDownLayers,'linear',2,'EndValues','nearest');
upInt=fillmissing(velUpLayers,'linear',2,'EndValues','nearest');
% For zenith data, do not use data that is too far away
smoothMins=60;
smoothPix=hcrTimeToPix(smoothMins,etime(datevec(data.time(2)),datevec(data.time(1))))+1;
maskZenithInt=movmean(velUpLayers,smoothPix,2,'omitnan');
upInt(isnan(maskZenithInt))=nan;

% Compare level altitudes between nadir and zenith and if they are close,
% calculate correction
zenithCorrectionRaw=[];
for ii=1:length(allMeanAltDown)
    for jj=1:length(allMeanAltUp)
        if abs(allMeanAltUp(jj)-allMeanAltDown(ii))<=3
            zenithCorrectionRaw=cat(1,zenithCorrectionRaw,downInt(ii,:)-upInt(jj,:));
        end
    end
end

% Average corrections over all levels
zenithCorrectionRaw=mean(zenithCorrectionRaw,1,'omitnan');
% If no correction was found, set it to zero everywhere
zenithCorrectionRaw(isnan(zenithCorrectionRaw))=0;
% Only up
zenithCorrectionRaw(data.elevation<=0)=nan;
% Smooth once more
zenithCorrection=movmean(zenithCorrectionRaw,smoothfactor,'omitnan');
zenithCorrection(isnan(zenithCorrectionRaw))=nan;
end