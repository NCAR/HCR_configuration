clear all
close all

infile='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/wband/save/20210529/20210529_191131.131_999_017.iwrf_ts';
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

    if feof(fileID)
        continue
    end

    if strcmp(packetID,'77770001')
        [syncIDthis,fileID]=readSyncIdData(fileID);
    elseif strcmp(packetID,'77770002')
        [radarInfo,fileID]=readRadarInfoData(fileID);
    elseif strcmp(packetID,'77770003')
        [scanSegment,fileID]=readScanSegmentData(fileID);
    elseif strcmp(packetID,'77770005')
        [tsProcessing,fileID]=readTsProcessingData(fileID);
    % elseif strcmp(packetID,'77770012')
    %     [statusXml,fileID]=readStatusXmlData(fileID);
    elseif strcmp(packetID,'77770008')
        [cal,fileID]=readCalData(fileID);
    elseif strcmp(packetID,'77770111')
        [geoRef,fileID]=readGeoRefData(fileID);
    elseif strcmp(packetID,'7777000C')
        [pulseHeader,IQ1,IQ2,range,fileID]=readPulseHeaderData(fileID);
    else
        dummy=fread(fileID,1,'uint8',lenBytes-56-1);
        warning(['Skipping packet ',packetID,'.'])
    end
end

fclose(fileID);