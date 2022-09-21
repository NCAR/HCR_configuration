function cloudID=makeCloudID(joinedFlag,minCloudSizePix)
% Create file with numbered clouds

cloudMask=joinedFlag==1;
cloudMask=bwareaopen(cloudMask,minCloudSizePix);

% Find connected
connMask=bwconncomp(cloudMask);

% Label
cloudID=double(labelmatrix(connMask));

end