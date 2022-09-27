function cloudID=makeCloudID(data,minCloudSizePix)
% Create file with numbered clouds
flagTemp=data.FLAG;

% Join over nadir/zenith edges
signElev=sign(data.elevation);
diffSign=diff(signElev);
switchInds=find(abs(diffSign)==2 | abs(diffSign)==1);

transMask=data.ANTFLAG==5;
transRegs=bwconncomp(transMask);

regList=[];
for ii=1:transRegs.NumObjects
    thisReg=transRegs.PixelIdxList{ii};
    isSwitch=intersect(thisReg,switchInds);
    if ~isempty(isSwitch) & flagTemp(18,thisReg(1)-1)==1 & flagTemp(18,thisReg(end)+1)==1
        flagTemp(18:end,thisReg)=nan;
        flagTemp(18,thisReg)=1;
        regList=cat(1,regList,thisReg);
    end
end

% Create number file
cloudMask=flagTemp==1;
cloudMask=bwareaopen(cloudMask,minCloudSizePix);

% Find connected
connMask=bwconncomp(cloudMask);

% Label
cloudID=double(labelmatrix(connMask));

% Fill in corner regions
for ii=1:length(regList)
    thisFlag=data.FLAG(:,regList(ii));
    firstNan=min(find(isnan(thisFlag)));
    if isempty(firstNan)
        firstNan=length(thisFlag);
    end
    cloudID(18:firstNan-1,regList(ii))=cloudID(18,regList(ii));
end
end