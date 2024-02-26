% Read HCR iwrf time series
% Documentation: https://github.com/NCAR/lrose-titan/blob/master/docs/pdf/IWRF_ts_format.pdf
% and https://github.com/NCAR/lrose-core/blob/master/codebase/libs/radar/src/include/radar/iwrf_data.h

clear all
close all

infile='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/wband/save/20210529/20210529_191131.131_999_017.iwrf_ts';
fileID=fopen(infile,'r','l');

guessSize=500000;
baseTime=datetime(1970,1,1);

sUTCall=nan(1,guessSize);
timeNSall=nan(1,guessSize);
data.latitude=nan(1,guessSize);
data.longitude=nan(1,guessSize);
data.altitude=nan(1,guessSize);
data.elevation=nan(1,guessSize);
data.azimuth=nan(1,guessSize);
data.pulse_width=nan(1,guessSize);
data.prt=nan(1,guessSize);

data.range=nan(770,guessSize);
data.IVc=nan(770,guessSize);
data.QVc=nan(770,guessSize);
data.IHx=nan(770,guessSize);
data.QHx=nan(770,guessSize);

beamCount=1;

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
    elseif strcmp(packetID,'77770012')
        % [statusXml,fileID]=readStatusXmlData(fileID);
        dummy=fread(fileID,1,'uint8',lenBytes-56-1);
    elseif strcmp(packetID,'77770008')
        [cal,fileID]=readCalData(fileID);
    elseif strcmp(packetID,'77770111')
        [geoRef,fileID]=readGeoRefData(fileID);
    elseif strcmp(packetID,'7777000C') % IQ data
        [pulseHeader,IQH,IQV,range,fileID]=readPulseHeaderData(fileID);
        % Add all vars
        sUTCall(:,beamCount)=sUTC;
        timeNSall(:,beamCount)=timeNS;
        data.range(:,beamCount)=range';
        data.latitude(:,beamCount)=geoRef.latitude;
        data.longitude(:,beamCount)=geoRef.longitude;
        data.altitude(:,beamCount)=geoRef.altMslKm/1000;
        data.elevation(:,beamCount)=pulseHeader.elevation;
        data.azimuth(:,beamCount)=pulseHeader.azimuth;
        data.pulse_width(:,beamCount)=tsProcessing.pulseWidthUS;
        data.prt(:,beamCount)=tsProcessing.prtUS;
        data.IVc(:,beamCount)=IQV(1,:)';
        data.QVc(:,beamCount)=IQV(2,:)';
        data.IHx(:,beamCount)=IQH(1,:)';
        data.QHx(:,beamCount)=IQH(2,:)';

        beamCount=beamCount+1;
    else
        dummy=fread(fileID,1,'uint8',lenBytes-56-1);
        warning(['Skipping packet ',packetID,'.'])
    end
end

fclose(fileID);

allNan=find(isnan(data.latitude));

allFields=fieldnames(data);
for ii=1:length(allFields)
    data.(allFields{ii})(:,allNan)=[];
end

sUTCall(:,allNan)=[];
timeNSall(:,allNan)=[];

data.time=baseTime+seconds(sUTCall+timeNSall*10^-9);

data.lambda=radarInfo.wavelengthCM/100;
data.dbz1km_v=cal.baseDbz1kmVC;
data.dbz1km_h=cal.baseDbz1kmHX;
data.noise_v=cal.noiseDbmVC;
data.noise_h=cal.noiseDbmHX;
data.rx_gain_v=cal.receiverGainDbVC;
data.rx_gain_h=cal.receiverGainDbHX;
data.beamwidth_v=cal.beamwidthDegV;
data.beamwidth_h=cal.beamwidthDegH;
