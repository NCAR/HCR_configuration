% Find stratiform and convective
function [stratConv maxAlt]= f_meltLayer_altOnly(data,meltLayer)

% Initialize output
stratConv=nan(size(data.time));

% Smooth data
dbzSmooth=movmean(data.dbzMasked,5,2);

% Find maximum reflectivity altitude
[dbzMax dbzMaxInd]=max(dbzSmooth,[],1);
maxInLin=sub2ind(size(data.DBZ),dbzMaxInd,1:length(data.time));
maxAlt=data.asl(maxInLin);

% Find melting layer altitude
meltInds=find(meltLayer>0);
meltAlt=data.asl(meltInds);

% Distance between max refl and melt alt
reflDist=maxAlt-meltAlt';
stratConv(reflDist>500)=1;
stratConv(reflDist<=500)=0;

% for ii=1:length(data.time)
%     dbzRay=dbzSmooth(:,ii);
%     altRay=data.asl(:,ii);
%     
% end
end