function [surfFlag cloudSum]=makeSurfFlag(data,linInd)
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

% Find cloud data
data.dbzMasked=data.DBZ;
data.dbzMasked(data.FLAG>1)=nan;

cloudSum=sum(~isnan(data.dbzMasked),1);
clearInds=cloudSum<3;
surfFlag(clearInds & isnan(surfFlag))=2;
    
cloudInds=any(~isnan(data.dbzMasked),1);
surfFlag(cloudInds & isnan(surfFlag))=1;
%surfFlag(find(any(~isnan(data.dbzMasked),1) & surfFlag~=0 & surfFlag~=2))=1;

% Remove missing surface echo
surfFlag(isnan(data.DBZ(linInd)))=0;

% Remove noise source cal, ant trans, and missing
surfFlag(find(any(data.FLAG>9,1) & surfFlag~=0))=0;

% Remove extinct
surfFlag(find(any(data.FLAG==3,1) & surfFlag~=0))=0;
end

