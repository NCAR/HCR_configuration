function [upRegs,upFrac,upMaxStrength,downMaxStrength,upMeanStrength,downMeanStrength]=upDownDraft(velIn,aslIn,rangePix,distance)
upRegs=[];

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
upMaxStrength=-(max(velSmooth(velSmooth<0),[],'omitnan'));
downMaxStrength=max(velSmooth(velSmooth>0),[],'omitnan');
upMeanStrength=-(mean(velSmooth(velSmooth<0),'omitnan'));
downMeanStrength=mean(velSmooth(velSmooth>0),'omitnan');

upFrac=sum(sum(velSmooth<0))/sum(sum(~isnan(velFilled)));

% close all
% surf(velSmooth,'EdgeColor','none');
% view(2)
% colormap('jet');
% caxis([-5,5]);

% Get updraft regions
updrafts=velSmooth<0;

updrafts=bwareaopen(updrafts,5);

upProps=regionprops('table',updrafts,'Area','BoundingBox','PixelIdxList');
upNum=size(upProps,1);

if upNum>0
    upRegWidth=upProps.BoundingBox(:,3).*distance/1000;
    upRegDepth=upProps.BoundingBox(:,4).*rangePix/1000;
    upRegPix=upProps.Area;
    upRegAsl=nan(upNum,1);
    upRegMean=nan(upNum,1);
    upRegMax=nan(upNum,1);
    for ii=1:upNum
        upRegAsl(ii)=mean(aslIn(upProps.PixelIdxList{ii}),'omitnan')/1000;
        upRegMean(ii)=-mean(velSmooth(upProps.PixelIdxList{ii}),'omitnan');
        upRegMax(ii)=-(min(velSmooth(upProps.PixelIdxList{ii}),[],'omitnan'));
    end
    upRegs=table(upProps.Area,upRegWidth,upRegDepth,upRegMean,upRegMax,upRegAsl,'VariableNames',{'numPix','width','depth','meanVel','maxVel','asl'});
end
end