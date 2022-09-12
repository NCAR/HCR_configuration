function cloudID=makeCloudID(joinedEcho,minCloudSizePix)
% Create file with numbered clouds

cloudMask=~isnan(joinedEcho);
cloudMask=bwareaopen(cloudMask,minCloudSizePix);

% Find connected
connMask=bwconncomp(cloudMask);

% Label
cloudID=double(labelmatrix(connMask));

end