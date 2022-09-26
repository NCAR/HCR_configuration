function cloudID=makeCloudID(data,minCloudSizePix)
% Create file with numbered clouds

% Join over nadir/zenith edges
signElev=sign(data.elevation);
diffSign=diff(signElev);
switchInds=find(abs(diffSign)==2 | abs(diffSign)==1);

transMask=data.ANTFLAG==5;











cloudMask=data.FLAG==1;
cloudMask=bwareaopen(cloudMask,minCloudSizePix);

% Find connected
connMask=bwconncomp(cloudMask);

% Label
cloudID=double(labelmatrix(connMask));

end