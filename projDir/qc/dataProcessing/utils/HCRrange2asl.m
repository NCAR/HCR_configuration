function asl = HCRrange2asl(range,elevation,altitude)
% Calculate m above see level from hcr range data
if size(elevation,1)~=1
    elevation=elevation';
    altitude=altitude';
end
if size(range,2)~=length(elevation)
    range=range';
end
asl=nan(size(range));
downInd=find(elevation<0);
upInd=find(elevation>=0);
asl(:,downInd)=-1*((range(:,downInd).*cosd(abs(elevation(downInd))-90))-altitude(downInd));
asl(:,upInd)=range(:,upInd).*cosd(abs(elevation(upInd))-90)+altitude(upInd);
end

