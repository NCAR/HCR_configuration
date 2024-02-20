function [pulseHeader,IQ1,IQ2,range,fileID]=readPulseHeaderData(fileID)

pulseHeader.seqNum=fread(fileID,1,'int64',4*4+2*4);

pulseHeader.elevation=fread(fileID,1,'float');
pulseHeader.azimuth=fread(fileID,1,'float');

pulseHeader.prt=fread(fileID,1,'float',4);
pulseHeader.pulseWidthUS=fread(fileID,1,'float');

pulseHeader.nGates=fread(fileID,1,'int32');

pulseHeader.nChannels=fread(fileID,1,'int32');
pulseHeader.iqEncoding=fread(fileID,1,'int32');
pulseHeader.hvFlag=fread(fileID,1,'int32');
pulseHeader.antTransition=fread(fileID,1,'int32',2*4);

pulseHeader.Ndata=fread(fileID,1,'int32');
pulseHeader.iqOffset=fread(fileID,[1,4],'int32');
pulseHeader.burstMag=fread(fileID,[1,4],'float');
pulseHeader.burstArg=fread(fileID,[1,4],'float');
pulseHeader.burstArgDiff=fread(fileID,[1,4],'float');

pulseHeader.scale=fread(fileID,1,'float');
pulseHeader.offset=fread(fileID,1,'float');

pulseHeader.NgatesBurst=fread(fileID,1,'int32');

pulseHeader.startRangeM=fread(fileID,1,'float');
pulseHeader.gateSpacingM=fread(fileID,1,'float');
pulseHeader.eventFlags=fread(fileID,1,'int32');
pulseHeader.txrxState=fread(fileID,1,'int32');
pulseHeader.rxPhaseDeg=fread(fileID,1,'float');

dummy=fread(fileID,[1,5],'int32');

%% IQ data
IQ1raw=fread(fileID,[2,pulseHeader.nGates],'int16');
IQ1=IQ1raw.*pulseHeader.scale+pulseHeader.offset;

IQ2raw=fread(fileID,[2,pulseHeader.nGates],'int16');
IQ2=IQ2raw.*pulseHeader.scale+pulseHeader.offset;

range=(0:pulseHeader.nGates-1).*pulseHeader.gateSpacingM+pulseHeader.startRangeM;
end