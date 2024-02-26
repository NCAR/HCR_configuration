function [radarInfo,fileID]=readRadarInfoData(fileID)

radarInfo.lat=fread(fileID,1,'float');
radarInfo.lon=fread(fileID,1,'float');
radarInfo.alt=fread(fileID,1,'float');

radarInfo.platform=fread(fileID,1,'int32');

radarInfo.beamwidthDegH=fread(fileID,1,'float');
radarInfo.beamwidthDegV=fread(fileID,1,'float');
radarInfo.wavelengthCM=fread(fileID,1,'float');

radarInfo.antGainH=fread(fileID,1,'float');
radarInfo.antGainV=fread(fileID,1,'float',25*4+64);

%dummy=fread(fileID,[1,25],'float');

%radarInfo.radarName=fread(fileID,[1,32],'*char');
%radarInfo.siteName=fread(fileID,[1,32],'*char');

end