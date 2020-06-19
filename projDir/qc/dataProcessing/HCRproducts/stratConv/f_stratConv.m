% Find stratiform and convective
function stratConv= f_meltLayer_altOnly(data,puzzle,meltLayer,meltArea)

disp('Finding stratiform and convective echo ...');

% Initialize output
stratConv=nan(size(data.DBZ));

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

% Loop through clouds
numClouds=max(max(puzzle));

for ii=1:numClouds
    stratConv1D=nan(size(data.time));
    
    cloudInd=find(puzzle==ii);
    [cloudR cloudC]=ind2sub(size(data.dbzMasked),cloudInd);
    cloudCu=unique(cloudC);
    
    % Max and min altitude of cloud    
%     cloudMask=nan(size(data.dbzMasked));
%     cloudMask(cloudInd)=data.dbzMasked(cloudInd);
    
    aslMask=nan(size(data.dbzMasked));
    aslMask(cloudInd)=data.asl(cloudInd);
    
    maxAltCloud=max(aslMask,[],1);
    minAltCloud=min(aslMask,[],1);
    
    % Shallow convection
    stratConv1D(meltAlt'>maxAltCloud)=2;
    
    % Stratiform because above melting layer
    stratConv1D(meltAlt'<minAltCloud)=0;   
    
    % Stratiform because melting layer detected
    meltYes=any(meltLayer==1,1);
    meltCloud=nan(size(data.time));
    meltCloud(cloudCu)=meltYes(cloudCu);
    stratConv1D(meltCloud==1)=0;
    
    % Calculate velocity steps
    % Detect steps in data
    velMask=nan(size(data.dbzMasked));
    velMask(cloudInd)=velSmooth(cloudInd);
    
    [velSteps,S1,S2] = ischange(velMask,1,'MaxNumChanges',1);
    
    [maxVel maxVelInd]=min(diff(velMask,1),[],1);
    
    % Find altitude of step
    stepInLin=find(velSteps==1);
    [stepInR,stepInC]=ind2sub(size(data.DBZ),stepInLin);
    
    stepAltLin=data.asl(stepInLin);
    
    stepAlt=nan(size(data.time));
    stepAlt(stepInC)=stepAltLin;
    
    % Distance between max vel and melt alt
    reflDist=stepAlt-meltAlt';
    
    % Check velocity data    
    stratConv1D(abs(reflDist)<100 & isnan(stratConv1D))=0;
    stratConv1D(reflDist>=100 & isnan(stratConv1D))=1;
    
    % Convective, no step in velocity found
    stratConv1D(isnan(stepAlt) & isnan(stratConv1D))=1;
    
    % Get rid of outliers
    stratConvMed=floor(movmedian(stratConv1D,9,'omitnan'));
    
    backMask=repmat(stratConvMed,size(aslMask,1),1);
    
    % Data above melt and below melt but not at melt
    noDataMelt=isnan(aslMask(meltInds));
    noMeltCols=find(noDataMelt' & meltAlt'>minAltCloud & meltAlt'<maxAltCloud);
    
%     noMeltMask=zeros(size(data.dbzMasked));
%     noMeltMask(:,noMeltCols)=1;
    
    for jj=1:length(noMeltCols)
        aslRay=aslMask(:,noMeltCols(jj));
        outRay=nan(size(aslRay));
        outRay(aslRay>meltAlt(noMeltCols(jj)))=0;
        outRay(aslRay<meltAlt(noMeltCols(jj)))=2;
        backMask(:,noMeltCols(jj))=outRay;
    end
            
    stratConv(~isnan(aslMask))=backMask(~isnan(aslMask));
end
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