% Analyze HCR clouds

clear all;
close all;

startTime=datetime(2019,8,25,16,30,0);
endTime=datetime(2019,8,25,19,0,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset
quality='field'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

figdir=['/h/eol/romatsch/hcrCalib/otherCalib/'];

directories.dataDir=HCRdir(project,quality,freqData);

%% Load data

data.eastward_velocity=[];
data.northward_velocity=[];
data.vertical_velocity=[];
data.pitch=[];
data.roll=[];
data.heading=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(directories.dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

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
%% Calculate

altFt=data.altitude* 3.28084;

groundSpd=sqrt(data.eastward_velocity.^2+data.northward_velocity.^2);
grSpdKt=groundSpd*1.94384;

vertFtS=data.vertical_velocity*3.28084;
%% Plot

close all

figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200]);

ax1=subplot(3,1,1);
plot(data.time,altFt./1000,'linewidth',2)
ylim([44 48]);
ylabel('Altitude (kFt)');

yyaxis right
plot(data.time,altFt./1000,'linewidth',2)
ylim([0 60]);
ylabel('Altitude (kFt)');
xlim([data.time(1),data.time(end)]);
grid on
title('Altitude')

ax2=subplot(3,1,2);
hold on
plot(data.time,grSpdKt,'linewidth',2)
ylim([400 520]);
ylabel('Ground speed (Knots)');
title('Ground speed and vertical velocity')

yyaxis right
plot(data.time,vertFtS,'linewidth',2);
ylim([-80 80]);
ylabel('Vert vel (Ft/s)');
xlim([data.time(1),data.time(end)]);
grid on

ax3=subplot(3,1,3);
hold on
plot(data.time,data.pitch*5,'linewidth',2)
plot(data.time,data.roll,'linewidth',2)
ylim([-40 40]);
ylabel('Pitch(*5) and roll angles (deg)');

yyaxis right
plot(data.time,data.heading,'linewidth',2)
ylim([0 400]);
ylabel('Heading (deg)');
xlim([data.time(1),data.time(end)]);
grid on

title('Pitch(*5), roll, and heading')

legend('Pitch','Roll','Heading');

set(gcf,'PaperPositionMode','auto')
print([figdir,'flightData_otrecRF09'],'-dpng','-r0');