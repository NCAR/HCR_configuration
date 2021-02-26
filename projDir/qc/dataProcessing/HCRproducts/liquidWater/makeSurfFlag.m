function surfFlag=makeSurfFlag(data,maxGate)
% Create flag of surface reflectivity
% 0 extinct or not usable
% 1 cloud
% 2 clear air 

surfFlag=nan(size(data.time));

%sort out non nadir pointing
surfFlag(data.elevation>-85)=0;

%sort out land
surfFlag(data.TOPO>0)=0;

% sort out data from below 2500m altitude
%surfFlag(data.altitude<2500)=0;

% Calculate reflectivity sum inside and outside ocean surface to
% distinguish clear air and cloud
reflTemp=data.DBZ;

% Remove bang
reflTemp(data.FLAG==6)=nan;

reflLin=10.^(reflTemp./10);
reflOceanLin=nan(size(data.time));
reflNoOceanLin=nan(size(data.time));

for ii=1:length(data.time)
    if (~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-5)) & ~isnan(maxGate(ii))
        reflRay=reflLin(:,ii);
        reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
        reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
    end
end

% Remove data where reflectivity outside of ocean swath is more than
% 0.8
clearAir=find(reflNoOceanLin<=0.5); % Is 0.8 in ssCal
surfFlag(clearAir)=2;
surfFlag(isnan(reflOceanLin))=0;

% Find cloud data
data.dbzMasked=data.DBZ;
data.dbzMasked(data.FLAG>1)=nan;
    
cloudInds=any(~isnan(data.dbzMasked),1);
surfFlag(cloudInds & isnan(surfFlag))=1;
%surfFlag(find(any(~isnan(data.dbzMasked),1) & surfFlag~=0 & surfFlag~=2))=1;

% Remove noise source cal, ant trans, and missing
surfFlag(find(any(data.FLAG>9,1) & surfFlag~=0))=0;

% Remove extinct
surfFlag(find(any(data.FLAG==3,1) & surfFlag~=0))=0;
end

