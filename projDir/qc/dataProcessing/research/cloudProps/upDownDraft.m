function [upNum,upFrac,upMaxWidth,upMaxDepth,upMaxStrength,downMaxStrength]=upDownDraft(velIn,rangePix,distance)
upMaxWidth=nan;
upMaxDepth=nan;

% Fill in gaps of missing data
allNan=all(isnan(velIn),1);
allNanMap=repmat(allNan,size(velIn,1),1);
nanSize=regionprops(allNan,'Area');
maxNan=max([nanSize(:).Area]);
if ~isempty(maxNan)
    velFilled=fillmissing(velIn,'linear',2,'EndValues','none','MaxGap',maxNan+1);
    velFilled(isnan(velIn) & allNanMap==0)=nan;
else
    velFilled=velIn;
end

% Smooth
velSmooth=smoothdata(velFilled,2,'movmedian',50,'omitnan');
velSmooth(isnan(velFilled))=nan;

% Overall stats
upMaxStrength=-(min(velSmooth(velSmooth<0),[],'omitnan'));
downMaxStrength=max(velSmooth(velSmooth>0),[],'omitnan');

upFrac=sum(sum(velSmooth<0))/sum(sum(~isnan(velFilled)));

% close all
% surf(velSmooth,'EdgeColor','none');
% view(2)
% colormap('jet');
% caxis([-5,5]);

% Get updraft regions
updrafts=velSmooth<0;

updrafts=bwareaopen(updrafts,5000);

upProps=regionprops('table',updrafts,'BoundingBox');
upNum=size(upProps,1);

if upNum>0
    upMaxWidth=max(upProps.BoundingBox(:,3))*distance/1000;
    upMaxDepth=max(upProps.BoundingBox(:,4))*rangePix/1000;
end
end