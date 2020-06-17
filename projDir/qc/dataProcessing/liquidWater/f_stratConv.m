% Find stratiform and convective
function [stratConv maxAlt]= f_meltLayer_altOnly(data,meltLayer,meltArea)

% Initialize output
stratConv=nan(size(data.time));

% Find melting layer altitude
meltInds=find(meltLayer>0);
meltAlt=data.asl(meltInds);

[meltInR meltInC]=ind2sub(size(data.DBZ),meltInds);

% Smooth data
%dbzSmooth=movmean(data.dbzMasked,10,2);

% Velocity
velMasked=data.VEL_CORR;
velMasked(data.FLAG>1)=nan;
velSmooth=movmean(velMasked,10,2);

velSmooth(:,data.elevation<0)=-velSmooth;

% Focus on melting layer
pixArea=round(meltArea/(data.range(2)-data.range(1)));

for ii=1:size(velSmooth,2)
    velSmooth(1:meltInR(ii)-pixArea,ii)=nan;
    velSmooth(meltInR(ii)+pixArea:end,ii)=nan;
end

[maxVel maxVelInd]=min(diff(velSmooth,1),[],1);

% for ii=1:size(velSmooth,2)
%     plot(velSmooth(:,ii),data.asl(:,ii));
% end

% Find maximum reflectivity altitude
maxInLin=sub2ind(size(data.DBZ),maxVelInd,1:length(data.time));
maxAlt=data.asl(maxInLin);

% Distance between max vel and melt alt
reflDist=maxAlt-meltAlt';
stratConv(abs(reflDist)<100)=0;
stratConv(reflDist>=100)=1;

% Is there data at melting layer?
meltData=data.dbzMasked(meltInds);

stratConv(isnan(meltData') & reflDist>0)=0;
stratConv(isnan(meltData') & reflDist<0)=2;

% % Find maximum reflectivity altitude
% [dbzMax dbzMaxInd]=max(dbzSmooth,[],1);
% maxInLin=sub2ind(size(data.DBZ),dbzMaxInd,1:length(data.time));
% maxAlt=data.asl(maxInLin);
% 
% % Find melting layer altitude
% meltInds=find(meltLayer>0);
% meltAlt=data.asl(meltInds);
% 
% % Distance between max refl and melt alt
% reflDist=maxAlt-meltAlt';
% stratConv(reflDist>500)=1;
% stratConv(reflDist<=500)=0;
% 
% % Is there data at melting layer?
% meltData=data.dbzMasked(meltInds);
% 
% stratConv(isnan(meltData') & reflDist>0)=0;
% stratConv(isnan(meltData') & reflDist<0)=2;

end