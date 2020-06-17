% Find stratiform and convective
function [stratConv stepAlt]= f_meltLayer_altOnly(data,meltLayer,meltArea)

% Initialize output
stratConv=nan(size(data.time));

% Find melting layer altitude
meltInds=find(meltLayer>0);
meltAlt=data.asl(meltInds);

[meltInR meltInC]=ind2sub(size(data.DBZ),meltInds);

% Mask velocity
velMasked=data.VEL_CORR;
velMasked(data.FLAG>1)=nan;

% Take care of down pointing
velMasked(:,data.elevation<0)=-velMasked;

% Focus on melting layer
pixArea=round(meltArea/(data.range(2)-data.range(1)));

for ii=1:size(velMasked,2)
    velMasked(1:meltInR(ii)-pixArea,ii)=nan;
    velMasked(meltInR(ii)+pixArea:end,ii)=nan;
end

% Smooth in range direction
velSmooth=movmedian(velMasked,10,1);

% Detect steps in data
[velSteps,S1,S2] = ischange(velSmooth,1,'MaxNumChanges',1);

% figure
% hold on
% for ii=1:size(velSmooth2,2)
%      plot(velSmooth2(:,ii),data.asl(:,ii));
%  end

[maxVel maxVelInd]=min(diff(velSmooth,1),[],1);

% Find altitude of step
stepInLin=find(velSteps==1);
[stepInR,stepInC]=ind2sub(size(data.DBZ),stepInLin);

stepAltLin=data.asl(stepInLin);

stepAlt=nan(size(data.time));
stepAlt(stepInC)=stepAltLin;

% Distance between max vel and melt alt
reflDist=stepAlt-meltAlt';
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