function [cal,fileID]=readCalData(fileID)

cal.wavelengthCM=fread(fileID,1,'float');
cal.beamwidthDegH=fread(fileID,1,'float');
cal.beamwidthDegV=fread(fileID,1,'float');

cal.antGainH=fread(fileID,1,'float');
cal.antGainV=fread(fileID,1,'float');

cal.pulseWidthUS=fread(fileID,1,'float');
cal.xmitPowerDbmH=fread(fileID,1,'float');
cal.xmitPowerDbmV=fread(fileID,1,'float');

cal.twoWayWaveguideLossDbH=fread(fileID,1,'float');
cal.twoWayWaveguideLossDbV=fread(fileID,1,'float');

cal.twoWayRadomeLossDbH=fread(fileID,1,'float');
cal.twoWayRadomeLossDbV=fread(fileID,1,'float');

cal.receiverMissmatchLossDb=fread(fileID,1,'float');

cal.radarConstantH=fread(fileID,1,'float');
cal.radarConstantV=fread(fileID,1,'float');

cal.noiseDbmHC=fread(fileID,1,'float');
cal.noiseDbmHX=fread(fileID,1,'float');
cal.noiseDbmVC=fread(fileID,1,'float');
cal.noiseDbmVx=fread(fileID,1,'float');

cal.receiverGainDbHC=fread(fileID,1,'float');
cal.receiverGainDbHX=fread(fileID,1,'float');
cal.receiverGainDbVC=fread(fileID,1,'float');
cal.receiverGainDbVX=fread(fileID,1,'float');

cal.baseDbz1kmHC=fread(fileID,1,'float');
cal.baseDbz1kmHX=fread(fileID,1,'float');
cal.baseDbz1kmVC=fread(fileID,1,'float');
cal.baseDbz1kmVX=fread(fileID,1,'float',4*4);

cal.noiseSourcePowerDbmH=fread(fileID,1,'float');
cal.noiseSourcePowerDbmV=fread(fileID,1,'float',10*4);

cal.receiverSlopeHC=fread(fileID,1,'float');
cal.receiverSlopeHX=fread(fileID,1,'float');
cal.receiverSlopeVC=fread(fileID,1,'float');
cal.receiverSlopeVX=fread(fileID,1,'float');

cal.i0DbmHC=fread(fileID,1,'float');
cal.i0DbmHX=fread(fileID,1,'float');
cal.i0DbmVC=fread(fileID,1,'float');
cal.i0DbmVX=fread(fileID,1,'float');

cal.dynamicRangeDbHC=fread(fileID,1,'float');
cal.dynamicRangeDbHX=fread(fileID,1,'float');
cal.dynamicRangeDbVC=fread(fileID,1,'float');
cal.dynamicRangeDbVX=fread(fileID,1,'float');

cal.kSquaredWater=fread(fileID,1,'float');
cal.dbzCorrection=fread(fileID,1,'float',49*4+32);

%dummy=fread(fileID,1,'float');

%radarInfo.radarName=fread(fileID,[1,32],'*char');
%radarInfo.siteName=fread(fileID,[1,32],'*char');

end