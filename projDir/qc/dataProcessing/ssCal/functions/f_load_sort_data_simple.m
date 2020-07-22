function [PLT, freq] = f_load_sort_data_simple(uniqueCases,indir)
%Function to determine the bias of HCR data and plot ocean scan data

% set up radar parameters, and partial radar constant term
% c = 3.0e8;
% tau = 2.56e-7;
% lambda = c/frq;
% w_dielec_sq = .69;
%
% temp = c * pi^5 * tau * w_dielec_sq /( 2*lambda^4 * 1e18);
% rc_terms = 10 * log10(temp);
%
% theta = 10.0;  % optimum angle for sigma0 determination
% dB_sigma10 = 6.0;   % expected sigma0 at 10 degrees

clear PLT;
PLT = struct();
PLT.time = [];
PLT.rota = [];
PLT.refl = [];  % max reflectivity
PLT.vpwr = [];
PLT.alt  = [];
PLT.vel  = [];
PLT.pitch = [];
PLT.elev = [];  % elevation angle is referenced to horiz plane;
% pitch and roll corrected.  -80 = 170 or 190 rotation
%PLT.angdif = []; % abs(rotation - 180) - (90 + elevation) (for downward)
PLT.roll = [];
PLT.range = [];
PLT.lon = [];
PLT.lat = [];
PLT.hdg = [];
PLT.pulseWidth=[];
reflInOrig=[];
xmitPowVin=[];
antGainVin=[];
beamWidthHin=[];
beamWidthVin=[];
waveGuideLossin=[];
radomeLossin=[];
recMismatchLossin=[];

freq=nan;

% make list with files to go through
startTime=datetime(uniqueCases(1:6));
endTime=datetime(uniqueCases(7:12));

oneCaseFileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
% got through all files and load data
for kk=1:size(oneCaseFileList,2);  % step through each input file in the given set
    inName=oneCaseFileList{kk};
    toLoc=strfind(inName,'to_');
    indataOrig=inName(1:toLoc);
    
    getInfile=dir([indataOrig '*']);
    
    if size(getInfile,1)>1
        disp('More than one file found. Taking first.');
    end
    indata=[getInfile(1).folder,'/',getInfile(1).name];
    
    % get frequency
    if kk==1
        freq=ncread(indata,'frequency');
    end
    
    % read in calibration data
    xmitPowVin=[xmitPowVin ncread(indata,'r_calib_xmit_power_v')];
    antGainVin=[antGainVin ncread(indata,'r_calib_antenna_gain_v')];
    beamWidthHin=[beamWidthHin ncread(indata,'radar_beam_width_h')];
    beamWidthVin=[beamWidthVin ncread(indata,'radar_beam_width_v')];
    waveGuideLossin=[waveGuideLossin ncread(indata,'r_calib_two_way_waveguide_loss_v')];
    radomeLossin=[radomeLossin ncread(indata,'r_calib_two_way_radome_loss_v')];
    recMismatchLossin=[recMismatchLossin ncread(indata,'r_calib_receiver_mismatch_loss')];
    
    % read in time and convert to datetime
    startTimeIn=ncread(indata,'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeRead=ncread(indata,'time')';
    
    timeIn=startTimeFile+seconds(timeRead);
    
    PLT.time = [ PLT.time;  timeIn' ];
    
    %read in data that won't change
    PLT.hdg   = [ PLT.hdg;   ncread(indata,'heading') ];  % hope for no heading near North
    PLT.lat   = [ PLT.lat;   ncread(indata,'latitude') ];
    PLT.lon   = [ PLT.lon;   ncread(indata,'longitude') ];
    PLT.pulseWidth   = [ PLT.pulseWidth;   ncread(indata,'pulse_width') ];
    
    %read in temporary data
    rotaIn=ncread(indata,'rotation');
    pitchIn=ncread(indata,'pitch');
    rollIn=ncread(indata,'roll');
    elevIn=ncread(indata,'elevation');
    rangeIn=ncread(indata,'range');
    velIn=ncread(indata,'VEL');
    reflIn=ncread(indata,'DBZ');
    pwrIn=ncread(indata,'DBMVC');
    altIn=ncread(indata,'altitude');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sort out bad data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reflInOrig=[reflInOrig,reflIn];
    
    % Low reflectivity values don't affect result
    reflIn(reflIn<-15)=nan;
    
    % Scans should only have reflectivity data at the very top and at the ocean surface.
    % We sort out scans with too many reflectivity values
    dataSize=sum(~isnan(reflIn),1);
    tooMuchInd=find(dataSize>50);
    
    reflIn(:,tooMuchInd)=nan;
    
    %sort out scans where roll angle is too large or small
    outRollInd=find(rollIn<-2 | rollIn>2);
    reflIn(:,outRollInd)=nan;
       
    % There are high reflectivities in the first gates that need to be removed.
    badRow=10;
    empty=false;
    
    % find first row that has less than 1% nans
    while (badRow<40 && ~empty)
        if length(find(~isnan(reflIn(badRow,:))))<size(reflIn,2)/100 | length(find(~isnan(reflIn(badRow,:))))<3;
            empty=true;
        end
        badRow=badRow+1;
    end;
    
    if badRow==40
        disp(['First row with only nans is above 40. Check data!']);
        reflIn(:,:)=nan;
    else
        reflIn(1:badRow-1,:)=nan;
    end
    
    [maxRefl maxGate]=nanmax(reflIn,[],1);
    
    %sort out scans with reflectivity values that are too high or low
    outRangeInd=find(maxRefl<0 | maxRefl>55);
    maxRefl(outRangeInd)=nan;
    
    elev=abs(elevIn+90);
%     %Double check elevation angle
%     elev = atand(sqrt((tand(pitchIn)).^2 + ...
%         (tand(rotaIn+ rollIn)./cosd(pitchIn)).^2));
%     
%     eldif = elevIn - elev + 90;
    
    %Get the linear index of the maximum reflectivity value
    maxGateLin=sub2ind(size(reflIn),maxGate,1:size(reflIn,2));
    
    %Check if ocean surface is at right altitude
    rangeMat=repmat(rangeIn,1,size(reflIn,2));
    sfcrng = rangeMat(maxGateLin);
    oceanSurf=sfcrng.*cosd(elev)';
    
    wrongAltInd=find(abs(altIn-oceanSurf')>100);
    maxRefl(wrongAltInd)=nan;
    
    % Power, range, and vel at maximum refl value
    maxpwr = pwrIn(maxGateLin);
    sfcvel = velIn(maxGateLin);
    
    %Get index where maximum reflectivity is empty and delete data
    emptyMax=find(isnan(maxRefl));
    if ~isempty(emptyMax)
        maxpwr(emptyMax)=nan;
        sfcrng(emptyMax)=nan;
        sfcvel(emptyMax)=nan;
    end
    
    PLT.refl = [ PLT.refl; maxRefl' ];
    PLT.vel  = [ PLT.vel;  sfcvel'];
    PLT.vpwr = [ PLT.vpwr; maxpwr' ];
    PLT.range = [ PLT.range; sfcrng' ];
    PLT.elev  = [ PLT.elev; elev ];
    %PLT.angdif = [ PLT.angdif; eldif ];
    PLT.rota = [ PLT.rota; rotaIn ];
    PLT.pitch = [ PLT.pitch; pitchIn ];
    PLT.roll = [ PLT.roll; rollIn ];
    PLT.alt = [ PLT.alt; altIn ];
end;

timeInds=find(PLT.time>startTime & PLT.time<endTime);

fieldsPLT=fields(PLT);

for jj=1:length(fieldsPLT)
    if ~isempty(PLT.(fieldsPLT{jj}))
        PLT.(fieldsPLT{jj})=PLT.(fieldsPLT{jj})(timeInds);
    end
end

reflInOrig=reflInOrig(:,timeInds);

%Warning when pitch angle is too large or small
outPitchInd=find(PLT.pitch<0 | PLT.pitch>4);
if ~isempty(outPitchInd)
    disp('Pitch angle is greater than 4 deg or smaller than 0 deg.');
end

xmitPowVu=unique(xmitPowVin);
if length(xmitPowVu)>1
    disp('xmitPowV has different values!!!');
end
PLT.xmitPowV=nanmean(xmitPowVu);

antGainVu=unique(antGainVin);
if length(antGainVu)>1
    disp('antGainV has different values!!!');
end
PLT.antGainV=antGainVu;

beamWidthHu=unique(beamWidthHin);
if length(beamWidthHu)>1
    disp('beamWidthH has different values!!!');
end
PLT.beamWidthH=beamWidthHu;

beamWidthVu=unique(beamWidthVin);
if length(beamWidthVu)>1
    disp('beamWidthV has different values!!!');
end
PLT.beamWidthV=beamWidthVu;

waveGuideLossu=unique(waveGuideLossin);
if length(waveGuideLossu)>1
    disp('waveGuideLoss has different values!!!');
end
PLT.waveGuideLoss=waveGuideLossu;

radomeLossu=unique(radomeLossin);
if length(radomeLossu)>1
    disp('radomeLoss has different values!!!');
end
PLT.radomeLoss=radomeLossu;

recMismatchLossu=unique(recMismatchLossin);
if length(recMismatchLossu)>1
    disp('recMismatchLoss has different values!!!');
end
PLT.recMismatchLoss=recMismatchLossu;

end

