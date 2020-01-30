function [PLT, freq] = f_load_sort_data(uniqueCases,indir)

clear PLT;
PLT = struct();
PLT.time = [];
PLT.rota = [];
PLT.refl = [];  % max reflectivity
PLT.vpwr = [];
PLT.hpwr = [];
PLT.alt  = [];
PLT.vel  = [];
PLT.pitch = [];
PLT.elev = [];  % elevation angle is referenced to horiz plane;
PLT.roll = [];
PLT.range = [];
PLT.rangeMat=[];
PLT.lon = [];
PLT.lat = [];
PLT.hdg = [];
PLT.pulseWidth=[];
PLT.dbz=[];
reflInOrig=[];
xmitPowVin=[];
antGainVin=[];
beamWidthHin=[];
beamWidthVin=[];
waveGuideLossin=[];
radomeLossin=[];
recMismatchLossin=[];
ksquaredin=[];
%wavelengthCmin=[];

freq=nan;

ID=[];

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
    
    % read in time and convert to datetime
    startTimeIn=ncread(indata,'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeRead=ncread(indata,'time')';
    timeIn=startTimeFile+seconds(timeRead);
    PLT.time = [ PLT.time, timeIn];
    
    rangeIn=ncread(indata,'range');
    rangeMat=repmat(rangeIn,1,length(timeIn));
    
    PLT.rota = [ PLT.rota, ncread(indata,'rotation')'];
    PLT.dbz = [PLT.dbz, ncread(indata,'DBZ')];
    PLT.vel = [PLT.vel, ncread(indata,'VEL')];
    PLT.vpwr = [PLT.vpwr, ncread(indata,'DBMVC')];
    PLT.hpwr = [PLT.hpwr, ncread(indata,'DBMHX')];
    PLT.alt  = [PLT.alt, ncread(indata,'altitude')'];
    PLT.pitch = [PLT.pitch, ncread(indata,'pitch')'];
    PLT.elev = [PLT.elev, ncread(indata,'elevation')'];  % elevation angle is referenced to horiz plane;
    PLT.roll = [PLT.roll, ncread(indata,'roll')'];
    PLT.range = [PLT.range, rangeMat];
    PLT.hdg   = [ PLT.hdg,   ncread(indata,'heading')' ];  % hope for no heading near North
    PLT.lat   = [ PLT.lat,   ncread(indata,'latitude')' ];
    PLT.lon   = [ PLT.lon,   ncread(indata,'longitude')' ];
    PLT.pulseWidth = [ PLT.pulseWidth, ncread(indata,'pulse_width')'];
    % read in calibration data
    xmitPowVin=[xmitPowVin ncread(indata,'r_calib_xmit_power_v')];
    antGainVin=[antGainVin ncread(indata,'r_calib_antenna_gain_v')];
    beamWidthHin=[beamWidthHin ncread(indata,'radar_beam_width_h')];
    beamWidthVin=[beamWidthVin ncread(indata,'radar_beam_width_v')];
    waveGuideLossin=[waveGuideLossin ncread(indata,'r_calib_two_way_waveguide_loss_v')];
    radomeLossin=[radomeLossin ncread(indata,'r_calib_two_way_radome_loss_v')];
    recMismatchLossin=[recMismatchLossin ncread(indata,'r_calib_receiver_mismatch_loss')];
    try
        ksquaredin=[ksquaredin ncread(indata,'r_calib_k_squared_water')];
    catch
        ksquaredin=[ksquaredin 0.711];
    end
    % wavelengthCmin=[wavelengthCmin varFromCfRadialString(indata,'wavelengthCm')];
    
    % If exist, read in model data
    
    if kk==1
        ncid = netcdf.open(indata,'nowrite');
        try
            ID = netcdf.inqVarID(ncid,'TEMP');
        end
        netcdf.close(ncid)
        if ~isempty(ID)
            PLT.topo=[];
            PLT.u=[];
            PLT.v=[];
            PLT.p=[];
            PLT.temp=[];
            PLT.rh=[];
            PLT.sst=[];
        end
    end
    if ~isempty(ID)
        PLT.topo = [PLT.topo, ncread(indata,'TOPO')'];
        PLT.u = [PLT.u, ncread(indata,'U_SURF')'];
        PLT.v = [PLT.v, ncread(indata,'V_SURF')'];
        PLT.p = [PLT.p, ncread(indata,'PRESS')];
        PLT.temp = [PLT.temp, ncread(indata,'TEMP')];
        PLT.rh = [PLT.rh, ncread(indata,'RH')];
        PLT.sst = [PLT.sst, ncread(indata,'SST')'];
    end
end

timeInds=find(PLT.time>startTime & PLT.time<endTime);

fieldsPLT=fields(PLT);

for jj=1:length(fieldsPLT)
    if ~isempty(PLT.(fieldsPLT{jj}))
        PLT.(fieldsPLT{jj})=PLT.(fieldsPLT{jj})(:,timeInds);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort out bad data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reflTemp=PLT.dbz;
PLT.reflMask=ones(size(PLT.time));
PLT.elev=abs(PLT.elev+90);

%sort out upward pointing
outElevInd=find(PLT.elev>90);
PLT.reflMask(outElevInd)=0;

% sort out data from below 2500m altitude
altInd=find(PLT.alt<2500);
PLT.reflMask(altInd)=0;

% Remove bang
reflTemp(1:17,:)=nan;

% Find ocean surface gate
[PLT.refl maxGate]=max(reflTemp,[],1);

% Exclude data with max gate=1
max1=find(maxGate==1);
PLT.reflMask(max1)=0;

%Get the linear index of the maximum reflectivity value
maxGateLin=sub2ind(size(reflTemp),maxGate,1:size(reflTemp,2));

%Check if ocean surface is at right altitude
sfcrng = PLT.range(maxGateLin);
oceanSurf=sfcrng.*cosd(PLT.elev);
wrongAltInd=find(abs(PLT.alt-oceanSurf)>100);
PLT.reflMask(wrongAltInd)=0;

% Calculate reflectivity sum inside and outside ocean surface
reflLin=10.^(reflTemp./10);
reflOceanLin=nan(size(PLT.time));
reflNoOceanLin=nan(size(PLT.time));
%PLT.reflNum=nan(size(PLT.time));

%reflTemp(reflTemp<-15)=nan;

for ii=1:length(PLT.time)
    if ~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-10)
        reflRay=reflLin(:,ii);
        reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
        reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
        %reflRayDBZ=reflTemp(:,ii);
        %PLT.reflNum(ii)=sum(~isnan(reflRayDBZ(1:maxGate(ii)-6)));
    end
end

% Remove data where reflectivity outside of ocean swath is more than
% 0.8
%reflTot=reflNoOceanLin+reflOceanLin;
%PLT.reflPercent=reflNoOceanLin./(reflTot).*100;
tooMuchRefl=find(reflNoOceanLin>0.8);
PLT.reflMask(tooMuchRefl)=0;
PLT.reflMask(isnan(reflOceanLin))=0;

% reflSumLin=sum(reflLin,1,'omitnan');
% PLT.reflSum=10*log10(reflSumLin);
% 
% % Remove data with too many reflectivity gates
% tooMany=find(PLT.reflNum>6);
% PLT.reflMask(tooMany)=0;
% % same for running mean
% tooManyMean=find(movmean(PLT.reflNum,50)>4);
% PLT.reflMask(tooManyMean)=0;

% remove data before and after drop outs
nanInds=find(isnan(PLT.refl));
if ~isempty(find(nanInds==1)) | ~isempty(find(nanInds==2)) | ~isempty(find(nanInds==3))
    nanInds=nanInds(4:end);
end
if ~isempty(find(nanInds==length(PLT.time))) | ~isempty(find(nanInds==length(PLT.time)-1)) | ~isempty(find(nanInds==length(PLT.time)-2))
    nanInds=nanInds(1:end-3);
end
PLT.reflMask(nanInds+1)=0;
PLT.reflMask(nanInds-1)=0;
PLT.reflMask(nanInds+2)=0;
PLT.reflMask(nanInds-2)=0;
PLT.reflMask(nanInds+3)=0;
PLT.reflMask(nanInds-3)=0;

% Power
PLT.maxpwrv = PLT.vpwr(maxGateLin);
PLT.maxpwrh = PLT.hpwr(maxGateLin);

PLT.vel = PLT.vel(maxGateLin);

PLT.refl(PLT.reflMask==0)=nan;

PLT.alt=PLT.alt';
PLT.pitch = PLT.pitch';
PLT.elev = PLT.elev';
PLT.roll = PLT.roll';
PLT.hdg   =  PLT.hdg';  % hope for no heading near North
PLT.lat   =  PLT.lat';
PLT.lon   =  PLT.lon';
PLT.pulseWidth = PLT.pulseWidth';
PLT.rota=PLT.rota';
PLT.rangeMat=PLT.range';
PLT.range=PLT.range(maxGateLin)';
PLT.range(wrongAltInd)=nan;
PLT.time=PLT.time';
PLT.reflMask=PLT.reflMask';
PLT.refl=PLT.refl';
%PLT.reflSum=PLT.reflSum';
%PLT.reflNum=PLT.reflNum';
%PLT.reflPercent=PLT.reflPercent';
PLT.reflOcean=reflOceanLin';
PLT.reflNoOcean=reflNoOceanLin';
PLT.maxpwrv=PLT.maxpwrv';
PLT.maxpwrh=PLT.maxpwrh';
PLT.vel=PLT.vel';

if ~isempty(ID)    
    PLT.topo = PLT.topo';
    PLT.u = PLT.u';
    PLT.v = PLT.v';
    PLT.p = PLT.p';
    PLT.temp = PLT.temp';
    PLT.rh = PLT.rh';
    PLT.sst = PLT.sst';
end

xmitPowVu=unique(xmitPowVin);
if length(xmitPowVu)>1
    disp('xmitPowV has different values!!!');
end
PLT.xmitPowV=mean(xmitPowVu,'omitnan');

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

ksquared=unique(ksquaredin);
if length(ksquared)>1
    disp('ksquared has different values!!!');
end
PLT.kSquaredWater=ksquared;

% wavelengthCm=unique(wavelengthCmin);
% if length(wavelengthCm)>1
%     disp('wavelength has different values!!!');
% end
% PLT.wavelengthCm=wavelengthCm;

end

