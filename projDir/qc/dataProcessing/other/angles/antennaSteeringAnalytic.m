% Actual antenna steering
clear all
close all

startTime=datetime(2019,8,7,16,0,0);
endTime=datetime(2019,8,7,16,30,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

dataDir=HCRdir(project,quality,freqData);

figdir='/h/eol/romatsch/hcrCalib/antennaSteering/analytic/';

%% Load data

data.tilt=[];
data.rotation=[];
data.roll=[];
data.pitch=[];
data.heading=[];
data.azimuth=[];
data.drift=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Plot pitch and roll

figure('Position',[200 500 1500 600],'DefaultAxesFontSize',12)

subplot(2,1,1)
hold on
plot(data.time,data.pitch)
legend('Pitch')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(2,1,2)
hold on
plot(data.time,data.roll)
legend('Roll')
xlim([data.time(1),data.time(end)])
ylabel('degree')

print([figdir,'pitchRoll_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
%% Calculate predicted tilt and rotation

tiltPred=nan(size(data.heading));
isZero=nan(size(data.heading));
rotationPred=nan(size(data.heading));

heading=data.heading;
roll=data.roll;
pitch=data.pitch;
azimuth=data.azimuth;
elevation=data.elevation;

x_t = sind(azimuth) .* cosd(elevation);
y_t = cosd(azimuth) .* cosd(elevation);
z_t = sind(elevation);

x_a = x_t .* (cosd(heading) .* cosd(roll) + sind(heading) .* sind(pitch) .* sind(roll)) + ...
    y_t .* (-sind(heading) .* cosd(roll) + cosd(heading) .* sind(pitch) .* sind(roll)) + ...
    z_t .* -cosd(pitch) .* sind(roll);

y_a = x_t .* sind(heading) .* cosd(pitch) + ...
    y_t .* cosd(heading) .* cosd(pitch) + ...
    z_t .* sind(pitch);

z_a = x_t .* (cosd(heading) .* sind(roll) - sind(heading) .* sind(pitch) .* cosd(roll)) + ...
    y_t .* (-sind(heading) .* sind(roll) - cosd(heading) .* sind(pitch) .* cosd(roll)) + ...
    z_t .* cosd(pitch) .* cosd(roll);

rotationPred=atan2d(x_a,z_a);
tiltPred=asind(y_a);

%% Plot

figure('Position',[200 500 1500 1000],'DefaultAxesFontSize',12)

subplot(4,1,1)
hold on
plot(data.time,data.rotation)
plot(data.time,rotationPred)
legend('RotationIn','RotationPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,2)
hold on
plot(data.time,data.tilt)
plot(data.time,tiltPred)
legend('TiltIn','TiltPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,3)
hold on
plot(data.time,data.elevation)
legend('ElevationIn')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,4)
hold on
plot(data.time,data.azimuth)
legend('AzimuthIn')
xlim([data.time(1),data.time(end)])
ylabel('degree')

print([figdir,'input_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')

%% Point nadir

wantedAz=nan(size(data.heading));
wantedAz(:)=0;
wantedElev=nan(size(data.heading));
wantedElev(data.elevation<0)=-90;
wantedElev(data.elevation>=0)=90;

tiltPred=nan(size(data.heading));
isZero=nan(size(data.heading));
rotationPred=nan(size(data.heading));

heading=data.heading;
roll=data.roll;
pitch=data.pitch;
azimuth=wantedAz;
elevation=wantedElev;

x_t = sind(azimuth) .* cosd(elevation);
y_t = cosd(azimuth) .* cosd(elevation);
z_t = sind(elevation);

x_a = x_t .* (cosd(heading) .* cosd(roll) + sind(heading) .* sind(pitch) .* sind(roll)) + ...
    y_t .* (-sind(heading) .* cosd(roll) + cosd(heading) .* sind(pitch) .* sind(roll)) + ...
    z_t .* -cosd(pitch) .* sind(roll);

y_a = x_t .* sind(heading) .* cosd(pitch) + ...
    y_t .* cosd(heading) .* cosd(pitch) + ...
    z_t .* sind(pitch);

z_a = x_t .* (cosd(heading) .* sind(roll) - sind(heading) .* sind(pitch) .* cosd(roll)) + ...
    y_t .* (-sind(heading) .* sind(roll) - cosd(heading) .* sind(pitch) .* cosd(roll)) + ...
    z_t .* cosd(pitch) .* cosd(roll);

rotationPred=atan2d(x_a,z_a);
tiltPred=asind(y_a);

%% Calculate elevation and azimuth

xsubt=cosd(rotationPred+data.roll) .* sind(data.drift) .* cosd(tiltPred) .* sind(data.pitch) + ...
    cosd(data.drift) .* sind(rotationPred+data.roll) .* cosd(tiltPred) - ...
    sind(data.drift) .* cosd(data.pitch) .* sind(tiltPred);

ysubt = -cosd(rotationPred+data.roll) .* cosd(data.drift) .* cosd(tiltPred) .* sind(data.pitch) + ...
    sind(data.drift) .* sind(rotationPred+data.roll) .* cosd(tiltPred) + ...
    cosd(data.pitch) .* cosd(data.drift) .* sind(tiltPred);

zsubt = cosd(data.pitch) .* cosd(tiltPred).* cosd(rotationPred+data.roll) + ...
    sind(data.pitch) .* sind(tiltPred);

xyz_field=[xsubt;ysubt;zsubt];
xyz_field=xyz_field;

az_t=atan2d(xyz_field(1,:),xyz_field(2,:));
elevReal=real(asind(xyz_field(3,:)));

track=data.heading+data.drift;
azReal=wrapTo360(az_t+track);

%% Plot

figure('Position',[200 500 1500 1000],'DefaultAxesFontSize',12)

subplot(4,1,1)
hold on
plot(data.time,data.rotation)
plot(data.time,rotationPred)
legend('RotationIn','RotationPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,2)
hold on
plot(data.time,data.tilt)
plot(data.time,tiltPred)
legend('TiltIn','TiltPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,3)
hold on
plot(data.time,data.elevation)
plot(data.time,elevReal)
legend('ElevationIn','ElevationPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,4)
hold on
plot(data.time,data.azimuth)
plot(data.time,azReal)
legend('AzimuthIn','AzimuthPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

print([figdir,'wanted_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')

%% Point somewhere else

% wantedAz=nan(size(data.heading));
% wantedAz(:)=70;

azMin=120;
azMax=170;
azStep=(azMax-azMin)/(length(data.heading)-1);
wantedAz=azMin:azStep:azMax;

elevMin=-90;
elevMax=90;
elevStep=(elevMax-elevMin)/(length(data.heading)-1);
wantedElev=elevMin:elevStep:elevMax;

tiltPred=nan(size(data.heading));
isZero=nan(size(data.heading));
rotationPred=nan(size(data.heading));

heading=data.heading;
roll=data.roll;
pitch=data.pitch;
azimuth=wantedAz;
elevation=wantedElev;

x_t = sind(azimuth) .* cosd(elevation);
y_t = cosd(azimuth) .* cosd(elevation);
z_t = sind(elevation);

x_a = x_t .* (cosd(heading) .* cosd(roll) + sind(heading) .* sind(pitch) .* sind(roll)) + ...
    y_t .* (-sind(heading) .* cosd(roll) + cosd(heading) .* sind(pitch) .* sind(roll)) + ...
    z_t .* -cosd(pitch) .* sind(roll);

y_a = x_t .* sind(heading) .* cosd(pitch) + ...
    y_t .* cosd(heading) .* cosd(pitch) + ...
    z_t .* sind(pitch);

z_a = x_t .* (cosd(heading) .* sind(roll) - sind(heading) .* sind(pitch) .* cosd(roll)) + ...
    y_t .* (-sind(heading) .* sind(roll) - cosd(heading) .* sind(pitch) .* cosd(roll)) + ...
    z_t .* cosd(pitch) .* cosd(roll);

rotationPred=atan2d(x_a,z_a);
tiltPred=asind(y_a);

%% Calculate elevation and azimuth

xsubt=cosd(rotationPred+data.roll) .* sind(data.drift) .* cosd(tiltPred) .* sind(data.pitch) + ...
    cosd(data.drift) .* sind(rotationPred+data.roll) .* cosd(tiltPred) - ...
    sind(data.drift) .* cosd(data.pitch) .* sind(tiltPred);

ysubt = -cosd(rotationPred+data.roll) .* cosd(data.drift) .* cosd(tiltPred) .* sind(data.pitch) + ...
    sind(data.drift) .* sind(rotationPred+data.roll) .* cosd(tiltPred) + ...
    cosd(data.pitch) .* cosd(data.drift) .* sind(tiltPred);

zsubt = cosd(data.pitch) .* cosd(tiltPred).* cosd(rotationPred+data.roll) + ...
    sind(data.pitch) .* sind(tiltPred);

xyz_field=[xsubt;ysubt;zsubt];
xyz_field=xyz_field;

az_t=atan2d(xyz_field(1,:),xyz_field(2,:));
elevReal=asind(xyz_field(3,:));

track=data.heading+data.drift;
azReal=wrapTo360(az_t+track);

%% Plot

figure('Position',[200 500 1500 1000],'DefaultAxesFontSize',12)

subplot(4,1,1)
hold on
plot(data.time,rotationPred)
legend('RotationPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,2)
hold on
plot(data.time,tiltPred)
legend('TiltPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,3)
hold on
plot(data.time,wantedElev)
plot(data.time,elevReal)
legend('ElevationWanted','ElevationPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

subplot(4,1,4)
hold on
plot(data.time,wantedAz)
plot(data.time,azReal)
legend('AzimuthWanted','AzimuthPred')
xlim([data.time(1),data.time(end)])
ylabel('degree')

print([figdir,'random_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')