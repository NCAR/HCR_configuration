% Calculate matrix output of Lee et al. 1994
% i.e. elevation and azimuth angles in earth coordinates

clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

drift_upper_limit=2;
drift_lower_limit=2;
tilt_error=0.05;
rotation_error=-0.25;

project='otrec'; % socrates, cset, aristo, otrec
quality='field'; % field, qc1, qc2
freq='10hz';

startTime=datetime(2019,7,31,18,1,0);
endTime=datetime(2019,7,31,18,11,0);

saveFig=1; % If 1, figures will be saved, if 0, figures will not be saved
outDir='/h/eol/romatsch/testOrder/';

% Desired variables. The variable name comies after the . and must be spelled exactly
% like in the CfRadial file

data.roll=[];
data.pitch=[];
data.heading=[];
data.rotation=[];
data.tilt=[];
data.drift=[];
data.northward_velocity=[];
data.eastward_velocity=[];
data.vertical_velocity=[];

data.VEL_RAW=[];
data.DBZ=[];

dataVars=fieldnames(data);

indir=HCRdir(project,quality,freq);

%% Load data
% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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

%% Process

% Convert to radian
drift_rad=deg2rad(data.drift);
pitch_rad=deg2rad(data.pitch);
roll_rad=deg2rad(data.roll);
theta_suba_rad=deg2rad(data.rotation-rotation_error);
tau_suba_rad=deg2rad(data.tilt-tilt_error);

% Compute y_t following equation 9 Lee et al (1994)
y_subt=-cos(theta_suba_rad+roll_rad).*cos(drift_rad).*cos(tau_suba_rad).*sin(pitch_rad)...
    +sin(drift_rad).*sin(theta_suba_rad+roll_rad).*cos(tau_suba_rad)...
    +cos(pitch_rad).*cos(drift_rad).*sin(tau_suba_rad);

% Compute z following equation 9 Lee et al (1994)
z=cos(pitch_rad).*cos(tau_suba_rad).*cos(theta_suba_rad+roll_rad)+sin(pitch_rad).*sin(tau_suba_rad);

% compute tau_t following equation 11 Lee et al (1994)
tau_subt=asin(y_subt);

% Compute phi following equation 17 Lee et al (1994)
phi=asin(z);

% Compute platform motion based on Eq 27 from Lee et al (1994)
ground_speed=sqrt(data.eastward_velocity.^2 + data.northward_velocity.^2);
vr_platform=-ground_speed.*sin(tau_subt)-data.vertical_velocity.*sin(phi);

VEL_corr=data.VEL_RAW-vr_platform;

% Compute the velocity of the ground
DBZ_tmp=data.DBZ;
DBZ_tmp(1:15,:)=0;
[bla ground_index]=max(DBZ_tmp,[],1);

for i=1:size(data.VEL_RAW,2)
    VEL_ground_raw(i)=data.VEL_RAW(ground_index(i),i);
    VEL_ground(i)=VEL_corr(ground_index(i),i);
end

figure;
plot(data.time,VEL_ground);
xlabel('Time, UTC');
ylabel('Corrected Ground Vr, m s^-^1');

figure;
plot(data.time,VEL_ground_raw,'r');
hold on;
plot(data.time,vr_platform,'k');
legend('Uncorrected Vr ground','Vr platform');
xlabel('Time, UTC');
ylabel('m s^-^1');

i_level=find(abs(data.roll<5) & data.elevation<-80 & abs(data.drift)<drift_upper_limit);
Vr_updown_mean=nanmean(VEL_ground(i_level));

if ~isempty(i_level)
    tilt_error_updown=asin(VEL_ground./ground_speed);
    tilt_error_deg_updown=tilt_error_updown*180/pi;
    tilt_error_updown_mean=nanmean(tilt_error_deg_updown(i_level));

    figure;
    hist(VEL_ground(i_level),50)
    xlabel('Ground Vr (m/s), upwind-downwind')
    title(strcat('Mean ground Vr = ',num2str(Vr_updown_mean),' Mean tilt error = ',num2str(tilt_error_updown_mean)))
end

    i_level2=find(abs(data.roll<5) & data.elevation<-80 & abs(data.drift)>drift_lower_limit);
    Vr_cross_mean=nanmean(VEL_ground(i_level2));
    
if ~isempty(i_level2)
    ratio=VEL_ground./(ground_speed.*sin(drift_rad));
    ratio(abs(ratio)>1)=nan;
    
    rotation_error_cross=-asind(ratio);
    rotation_error_cross_mean=nanmean(rotation_error_cross(i_level2));
    
    figure;
    hist(VEL_ground(i_level2),50)
    xlabel('Ground Vr (m/s), cross wind')
    title(strcat('Mean ground Vr = ',num2str(Vr_cross_mean),' Mean rotation error = ',num2str(rotation_error_cross_mean)))   
end

