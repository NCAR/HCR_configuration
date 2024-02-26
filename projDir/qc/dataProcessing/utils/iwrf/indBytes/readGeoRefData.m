function [geoRef,fileID]=readGeoRefData(fileID)

geoRef.unitNum=fread(fileID,1,'int32',4);

geoRef.altMslKm=fread(fileID,1,'float',4);

geoRef.EWvelMpS=fread(fileID,1,'float');
geoRef.NSvelMpS=fread(fileID,1,'float');
geoRef.VertVelMpS=fread(fileID,1,'float');

geoRef.headingDeg=fread(fileID,1,'float');
geoRef.rollDeg=fread(fileID,1,'float');
geoRef.pitchDeg=fread(fileID,1,'float');
geoRef.driftDeg=fread(fileID,1,'float');
geoRef.rotationDeg=fread(fileID,1,'float');
geoRef.tiltDeg=fread(fileID,1,'float',7*4);

geoRef.longitude=fread(fileID,1,'float64');
geoRef.latitude=fread(fileID,1,'float64',2*4+24*4);

end