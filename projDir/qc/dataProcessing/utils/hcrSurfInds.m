function [linInd rowInd rangeToSurf] = hcrSurfInds(data)
rangeTemp=data.range;
DBZmask=data.DBZ; %Mask with DBZ values that are close to the surface

if isfield(data,'topo')
    compAlt=data.altitude-data.topo;
    oceanInd=zeros(size(compAlt));
    oceanInd(data.topo<1)=1;
elseif isfield(data,'TOPO')
    compAlt=data.altitude-data.TOPO;
    oceanInd=zeros(size(compAlt));
    oceanInd(data.TOPO<1)=1;
else
    compAlt=data.altitude;
    oceanInd=zeros(size(compAlt));
end

rangeTemp(find(abs(data.range.*cosd(data.elevation+90)-compAlt)>100 & oceanInd==1))=nan;
rangeTemp(find(abs(data.range.*cosd(data.elevation+90)-compAlt)>500 & oceanInd==0))=nan;

DBZmask(isnan(rangeTemp))=nan;

upInd=find(data.elevation>-45);
DBZmask(:,upInd)=nan;

[bla ground_index]=max(DBZmask,[],1);
wrong_ground_ind=find(ground_index==1);
rowInd=ground_index;
rowInd(wrong_ground_ind)=nan;

% Mask out when we don't see the whole ground swath
indOut=zeros(size(oceanInd));
for ii=1:length(rowInd)
    if ~isnan(rowInd(ii)) & compAlt(ii)>100
        rayDBZ=data.DBZ(rowInd(ii)-5:min(rowInd(ii)+5,size(DBZmask,1)),ii);
        if length(find(~isnan(rayDBZ)))<5
            indOut(ii)=1;
        end
    end
end

% Convert to linear indices
linInd=sub2ind(size(data.DBZ),ground_index,1:length(data.altitude));

rangeToSurf=data.range(linInd); %Range to surface, should be very close to alt
tooLow=find(rangeToSurf<100);

indOut(tooLow)=1;

% Surf ind shouldn't be too close to end of range
indOut(rowInd>size(DBZmask,1)-5)=1;

rangeToSurf(indOut==1)=nan;
rowInd(indOut==1)=nan;
end

