function [upRegs,upFrac,upMaxStrength,downMaxStrength,upMeanStrength,downMeanStrength]=upDownDraft(velIn,aslIn,rangePix,distance,lon,lat)
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

% % Smooth
% velSmooth=smoothdata(velFilled,2,'movmedian',50,'omitnan');
% velSmooth(isnan(velFilled))=nan;

% Overall stats
upMaxStrength=-(min(velFilled(velFilled<0),[],'omitnan'));
if isempty(upMaxStrength)
    upMaxStrength=nan;
end
downMaxStrength=max(velFilled(velFilled>0),[],'omitnan');
if isempty(downMaxStrength)
    downMaxStrength=nan;
end
upMeanStrength=-(mean(velFilled(velFilled<0),'omitnan'));
downMeanStrength=mean(velFilled(velFilled>0),'omitnan');

upFrac=sum(sum(velFilled<0))/sum(sum(~isnan(velFilled)));

% Get updraft regions
updrafts=velFilled<0;

updrafts=bwareaopen(updrafts,5);

upProps=regionprops('table',updrafts,'Area','BoundingBox','PixelIdxList','Centroid');
upNum=size(upProps,1);

if upNum>0
    upRegWidth=upProps.BoundingBox(:,3).*distance/1000;
    upRegDepth=upProps.BoundingBox(:,4).*rangePix/1000;
    pixKM2=distance/1000*rangePix/1000;
    upRegArea=upProps.Area*pixKM2;
    upRegCloudAltPerc=nan(upNum,1);
    upRegMean=nan(upNum,1);
    upRegMax=nan(upNum,1);
    upRegLon=lon(round(upProps.Centroid(:,1)))';
    upRegLat=lat(round(upProps.Centroid(:,1)))';

    % Normalized cloud altitude
    aslScaled=aslIn-max(aslIn(:));
    aslNorm=aslScaled./min(aslScaled(:));

    for ii=1:upNum
        upRegCloudAltPerc(ii)=mean(aslNorm(upProps.PixelIdxList{ii}),'omitnan')*100;
        upRegMean(ii)=-mean(velFilled(upProps.PixelIdxList{ii}),'omitnan');
        upRegMax(ii)=-(min(velFilled(upProps.PixelIdxList{ii}),[],'omitnan'));
    end
    upRegs=table(upRegArea,upRegWidth,upRegDepth,upRegMean,upRegMax,upRegCloudAltPerc,upRegLon,upRegLat, ...
        'VariableNames',{'area','width','depth','meanVel','maxVel','cloudAltPerc','lon','lat'});
end
end