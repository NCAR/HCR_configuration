function [tsProcessing,fileID]=readTsProcessingData(fileID)
tsProcessing.dummy=fread(fileID,1,'int32',3*4);
tsProcessing.prtUS=fread(fileID,1,'float');
tsProcessing.prtUS2=fread(fileID,1,'float',4+4);
tsProcessing.pulseWidthUS=fread(fileID,1,'float');
tsProcessing.startRange=fread(fileID,1,'float');
tsProcessing.gateSpacing=fread(fileID,1,'float',4*4+3*4);
tsProcessing.polMode=fread(fileID,1,'int32',+4*4+2*4);
tsProcessing.numPrt=fread(fileID,1,'int32');
tsProcessing.prtUS3=fread(fileID,1,'float');
tsProcessing.prtUS4=fread(fileID,1,'float');
tsProcessing.blockModePrt2=fread(fileID,1,'int32');
tsProcessing.blockModePrt3=fread(fileID,1,'int32');
tsProcessing.blockModePrt4=fread(fileID,1,'int32');
tsProcessing.polSyncMode=fread(fileID,1,'uint32',4*18);
end