function [scanSegment,fileID]=readScanSegmentData(fileID)
scanSegment.scanMode=fread(fileID,1,'int32',4*4+11*4+3*4+512*4+5*4+3*4+455*4+32*2);
end