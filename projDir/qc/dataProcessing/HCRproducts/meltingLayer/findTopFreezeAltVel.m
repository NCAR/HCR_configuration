function [zeroAlt,zeroVel]=findTopFreezeAltVel(data)
% Find all freezing levels
oddAngles=find(data.elevation>-70 & data.elevation<70);
upInds=find(data.elevation>=0);
downInds=find(data.elevation<0);

data.TEMP(:,oddAngles)=nan;
data.TEMP(data.TEMP==0)=-0.01;
signTemp=sign(data.TEMP);
zeroDeg=diff(signTemp,1);

zeroDeg(:,upInds)=-zeroDeg(:,upInds);

zeroDeg=cat(1,zeroDeg,nan(size(data.time)));
zeroDeg(zeroDeg>0 & ~isnan(zeroDeg))=1;

allNan=all(isnan(zeroDeg),1);
findOne=any(zeroDeg==1,1);

belowSurfInds=find(allNan==0 & findOne==0);

indMat=repmat((1:size(data.DBZ_MASKED,1))',1,size(data.DBZ_MASKED,2));
indMat(zeroDeg~=1)=nan;

% Find highest freezing levels
% Zenith
indUp=indMat(:,upInds);
zeroIndUp=max(indUp,[],1,'omitnan');
% Nadir
indDown=indMat(:,downInds);
zeroIndDown=min(indDown,[],1,'omitnan');

zeroInd=nan(size(data.time));
zeroInd(upInds)=zeroIndUp;
zeroInd(downInds)=zeroIndDown;
zeroIndNoNan=zeroInd;
colZero=1:length(zeroInd);
colZero(isnan(zeroInd))=[];
zeroIndNoNan(isnan(zeroInd))=[];
linZeroIndNoNan=sub2ind(size(data.DBZ_MASKED),zeroIndNoNan,colZero);

% Altitude
zeroAltNonNan=data.asl(linZeroIndNoNan);
zeroAlt=nan(size(data.time));
zeroAlt(colZero)=zeroAltNonNan;
zeroAlt(belowSurfInds)=data.TOPO(belowSurfInds);

% Velocity
zeroVelNonNan=data.VEL_MASKED(linZeroIndNoNan);
zeroVel=nan(size(data.time));
zeroVel(colZero)=zeroVelNonNan;
zeroVel(upInds)=-zeroVel(upInds);
end