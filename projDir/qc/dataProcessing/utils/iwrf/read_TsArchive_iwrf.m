clear all
close all

infile='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/wband/save/20210529/20210529_153908.125_999_060.iwrf_ts';
fileID=fopen(infile,'r','l');

% Read data
while ~feof(fileID) % Run until the end of file
    
    packetIDin=fread(fileID,1,'int32');
    lenBytes=fread(fileID,1,'int32');
    seqNum=fread(fileID,1,'int64');
    versionNum=fread(fileID,1,'int32');
    radarID=fread(fileID,1,'int32');
    sUTC=fread(fileID,1,'int64');
    timeNS=fread(fileID,1,'int32');
    res=fread(fileID,[1,5],'int32');

    packetID=dec2hex(packetIDin);

    if strcmp(packetID,'77770001')
        [syncIDthis,fileID]=readSyncIdData(fileID);
    elseif strcmp(packetID,'77770002')
        [radarInfo,fileID]=readRadarInfoData(fileID);
    elseif strcmp(packetID,'77770003')
        [scanSegment,fileID]=readScanSegmentData(fileID);
    elseif strcmp(packetID,'77770005')
        [tsProcessing,fileID]=readTsProcessingData(fileID);
    elseif strcmp(packetID,'77770012')
        [statusXml,fileID]=readStatusXmlData(fileID);
    end
end

fclose(fileID);