% Read HCR iwrf time series
% Documentation: https://github.com/NCAR/lrose-titan/blob/master/docs/pdf/IWRF_ts_format.pdf
% and https://github.com/NCAR/lrose-core/blob/master/codebase/libs/radar/src/include/radar/iwrf_data.h

clear all
close all
tic
infile='/scr/virga1/rsfdata/projects/spicule/hcr/time_series/wband/save/20210529/20210529_191131.131_999_017.iwrf_ts';
fileID=fopen(infile,'r','l');

% Read data
disp(['Reading file ',infile]);
dataAll=fread(fileID, [1,inf], '*uint8');
fclose(fileID);

numGates=770;
guessSize=ceil(length(dataAll)/numGates/8);
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

data.range=nan(numGates,guessSize);
data.IVc=nan(numGates,guessSize);
data.QVc=nan(numGates,guessSize);
data.IHx=nan(numGates,guessSize);
data.QHx=nan(numGates,guessSize);

beamCount=1;

%% Organize data
disp('Organizing data ...')
ii=1;
while ii<length(dataAll) % Run until the end of file
    packetIDin=typecast(dataAll(ii:ii+3),'int32');
    ii=ii+4;
    lenBytes=typecast(dataAll(ii:ii+3),'int32');
    ii=ii+20;
    sUTC=typecast(dataAll(ii:ii+7),'int64');
    ii=ii+8;
    timeNS=typecast(dataAll(ii:ii+3),'int32');
    ii=ii+24;
  
    packetID=dec2hex(packetIDin);

    if strcmp(packetID,'77770001')
        ii=ii+lenBytes-56;
    elseif strcmp(packetID,'77770002')
        ii=ii+24;
        wavelengthCM=typecast(dataAll(ii:ii+3),'single');
        ii=ii+176;
    elseif strcmp(packetID,'77770003')
        ii=ii+lenBytes-56;
    elseif strcmp(packetID,'77770005')
        ii=ii+16;
        prtUS=typecast(dataAll(ii:ii+3),'single');
        ii=ii+12;
        pulseWidthUS=typecast(dataAll(ii:ii+3),'single');
        ii=ii+172;
    elseif strcmp(packetID,'77770012')
        ii=ii+lenBytes-56;
    elseif strcmp(packetID,'77770008')
        ii=ii+4;
        beamwidthDegH=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        beamwidthDegV=typecast(dataAll(ii:ii+3),'single');
        ii=ii+56;
        noiseDbmHX=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        noiseDbmVC=typecast(dataAll(ii:ii+3),'single');
        ii=ii+12;
        receiverGainDbHX=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        receiverGainDbVC=typecast(dataAll(ii:ii+3),'single');
        ii=ii+12;
        baseDbz1kmHX=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        baseDbz1kmVC=typecast(dataAll(ii:ii+3),'single');
        ii=ii+356;
    elseif strcmp(packetID,'77770111')
        ii=ii+8;
        altMslKm=typecast(dataAll(ii:ii+3),'single');
        ii=ii+72;
        longitude=typecast(dataAll(ii:ii+7),'double');
        ii=ii+8;
        latitude=typecast(dataAll(ii:ii+7),'double');
        ii=ii+112;
    elseif strcmp(packetID,'7777000C') % IQ data
        ii=ii+32;
        elevation=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        azimuth=typecast(dataAll(ii:ii+3),'single');
        ii=ii+16;
        nGates=typecast(dataAll(ii:ii+3),'int32');
        ii=ii+96;
        scale=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        offset=typecast(dataAll(ii:ii+3),'single');
        ii=ii+8;
        startRangeM=typecast(dataAll(ii:ii+3),'single');
        ii=ii+4;
        gateSpacingM=typecast(dataAll(ii:ii+3),'single');
        ii=ii+36;

        % Create I/Q
        IQHraw=nan(2,nGates);
        IQHraw(1,:)=typecast(dataAll(ii:ii+2*nGates-1),'int16');
        ii=ii+2*nGates;
        IQHraw(2,:)=typecast(dataAll(ii:ii+2*nGates-1),'int16');
        ii=ii+2*nGates;
        IQH=IQHraw.*scale+offset;

        IQVraw=nan(2,nGates);
        IQVraw(1,:)=typecast(dataAll(ii:ii+2*nGates-1),'int16');
        ii=ii+2*nGates;
        IQVraw(2,:)=typecast(dataAll(ii:ii+2*nGates-1),'int16');
        ii=ii+2*nGates;
        IQV=IQVraw.*scale+offset;

        range=single(0:nGates-1).*gateSpacingM+startRangeM+gateSpacingM/2;

        % Add all vars
        sUTCall(:,beamCount)=sUTC;
        timeNSall(:,beamCount)=timeNS;
        data.range(:,beamCount)=range';
        data.latitude(:,beamCount)=latitude;
        data.longitude(:,beamCount)=longitude;
        data.altitude(:,beamCount)=altMslKm/1000;
        data.elevation(:,beamCount)=elevation;
        data.azimuth(:,beamCount)=azimuth;
        data.pulse_width(:,beamCount)=pulseWidthUS;
        data.prt(:,beamCount)=prtUS;
        data.IVc(:,beamCount)=IQV(1,:)';
        data.QVc(:,beamCount)=IQV(2,:)';
        data.IHx(:,beamCount)=IQH(1,:)';
        data.QHx(:,beamCount)=IQH(2,:)';

        beamCount=beamCount+1;
    else
        ii=ii+lenBytes-56;
        warning(['Skipping packet ',packetID,'.'])
    end
end

allNan=find(isnan(data.latitude));

allFields=fieldnames(data);
for ii=1:length(allFields)
    data.(allFields{ii})(:,allNan)=[];
end

sUTCall(:,allNan)=[];
timeNSall(:,allNan)=[];

data.time=baseTime+seconds(sUTCall+timeNSall*10^-9);

data.lambda=wavelengthCM/100;
data.dbz1km_v=baseDbz1kmVC;
data.dbz1km_h=baseDbz1kmHX;
data.noise_v=noiseDbmVC;
data.noise_h=noiseDbmHX;
data.rx_gain_v=receiverGainDbVC;
data.rx_gain_h=receiverGainDbHX;
data.beamwidth_v=beamwidthDegV;
data.beamwidth_h=beamwidthDegH;
toc