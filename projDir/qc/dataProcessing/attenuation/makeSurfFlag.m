function [surfFlag,cloudSum]=makeSurfFlag(data,linInd)
% Create flag of surface reflectivity
% 0 extinct
% 1 cloud
% 2 clear air 

surfFlag=nan(size(data.time));

%sort out land
surfFlag(data.TOPO>0)=nan;

% Find cloud data
cloudSum=sum(~isnan(data.dbzMasked),1);
surfFlag(cloudSum>0)=1;

% Find clear data
surfFlag(~isnan(data.DBZ(linInd)) & isnan(surfFlag))=2;

% 
% Remove extinct
%surfFlag(find(any(data.FLAG==3,1)))=nan;
end

