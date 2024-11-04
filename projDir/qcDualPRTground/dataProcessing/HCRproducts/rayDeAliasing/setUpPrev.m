function [velPrev,prevCount,prevKeep,flipYes]=setUpPrev(finalRay,velPrev,prevCount,prevKeep,flipYes,elev,dataFreq,defaultPrev)
% Set up previous ray for time consistency check

% Check if pointing direction is changed
if abs(elev)<20
    flipYes=1;
end
if flipYes
    velPrev=repmat(defaultPrev,length(finalRay),1);
    if abs(elev)>89
        flipYes=0;
    end
else
    % Check if elevation angle is wrong
    if abs(elev)<60
        velForPrevMask=logical(zeros(size(velPrev)));
        velPrev=nan(size(velPrev));
    else
        % Decide which velocities to keep for time consistency check: join
        % close regions and remove small isolated regions
        velForPrev=finalRay;
        velForPrev=movmean(velForPrev,5,'omitnan');
        velForPrev=movmean(velForPrev,5,'includenan');
        velForPrevMask=~isnan(velForPrev);
        velForPrevMask=bwareaopen(velForPrevMask,15);
        velForPrevMask(~isnan(velPrev))=1;

        % Create new velPrev
        velPrev=movmedian(finalRay,35,'omitnan'); % Was 20
        velPrev(isnan(finalRay))=nan;
        velPrev(~velForPrevMask)=nan;
    end

    % Add new values to prevKeep and handle counts
    prevKeep(velForPrevMask)=velPrev(velForPrevMask);
    prevCount(velForPrevMask)=0;
    prevCount(~velForPrevMask)=prevCount(~velForPrevMask)+1;
    prevKeep(prevCount>=dataFreq*5)=nan;

    % Add old velocities
    velPrev(isnan(velPrev))=prevKeep(isnan(velPrev));

    % Interpolate prev
    prevMedLarge=movmedian(velPrev,75,'omitnan'); % Was 50

    velPrev(isnan(prevMedLarge))=defaultPrev;
    velPrev=fillmissing(velPrev,'linear','EndValues','nearest');
end
end