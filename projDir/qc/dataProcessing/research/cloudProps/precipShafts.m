function precShafts=precipShafts(shaftMap,dbzMap,velMap,aslMap,groundDist)

% Fill in gaps of missing data
allNan=all(isnan(dbzMap),1);
allNanMap=repmat(allNan,size(dbzMap,1),1);
nanSize=regionprops(allNan,'Area');
maxNan=max([nanSize(:).Area]);
if ~isempty(maxNan)
    dbzFilled=fillmissing(dbzMap,'linear',2,'EndValues','none','MaxGap',maxNan+1);
    dbzFilled(isnan(dbzMap) & allNanMap==0)=nan;
else
    dbzFilled=dbzMap;
end

% Fill in gaps of missing data
allNan=all(isnan(velMap),1);
allNanMap=repmat(allNan,size(velMap,1),1);
nanSize=regionprops(allNan,'Area');
maxNan=max([nanSize(:).Area]);
if ~isempty(maxNan)
    velFilled=fillmissing(velMap,'linear',2,'EndValues','none','MaxGap',maxNan+1);
    velFilled(isnan(velMap) & allNanMap==0)=nan;
else
    velFilled=velMap;
end

% Smooth
velSmooth=smoothdata(velFilled,2,'movmedian',50,'omitnan');
velSmooth(isnan(velFilled))=nan;

% Get shaft
minAsl=min(aslMap(:),[],'omitnan');

shaftMap(aslMap>minAsl+200)=0;

shaft1D=sum(shaftMap,1);
shaft1D=shaft1D>0;

shaftKM=sum(shaft1D)*groundDist/1000;
shaftFrac=sum(shaft1D)/size(shaftMap,2);
shaftMeanRefl=mean(dbzFilled(shaftMap==1));
shaftMaxRefl=max(dbzFilled(shaftMap==1),[],'omitnan');
shaftMeanVel=mean(velSmooth(shaftMap==1));
shaftMaxVel=max(velSmooth(shaftMap==1),[],'omitnan');

precShafts=table(shaftKM,shaftFrac,shaftMeanRefl,shaftMaxRefl,shaftMeanVel,shaftMaxVel, ...
    'VariableNames',{'shaftKM','frac','meanRef','maxRefl','meanVel','maxVel'});
end