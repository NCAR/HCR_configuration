% Calculate PID from HCR HSRL combined data

clear all
%close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='2hzMerged'; % 10hz, 100hz, or 2hz

startTime=datetime(2018,1,29,01,40,0); %Wang_Rauber
endTime=datetime(2018,1,29,01,59,0); %Wang_Rauber

ylimits=[0 1.5];

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/pid_hcr/',project,'Mat/'];

infilePID=[indir,'era5.pid.20180128_224515_to_20180129_060956.Flight6.mat'];
infileTime=[indir,'era5.time.20180128_224515_to_20180129_060956.Flight6.mat'];
%infilePID=[indir,'era5.pid.20180129_013000_to_20180129_020000.Flight6.mat'];
%infileTime=[indir,'era5.time.20180129_013000_to_20180129_020000.Flight6.mat'];

%HCR data
PIDin=load(infilePID);
PID=PIDin.pid;
Timein=load(infileTime);
time=Timein.timeHCR;

gates=0:20:(size(PID,1)-1)*20;
asl=repmat(gates',1,size(PID,2));
asl=asl./1000;

timeInds=find(time>=startTime & time<=endTime);

time=time(timeInds);
asl=asl(:,timeInds);
PID=PID(:,timeInds);

cscale_comb=[1,1,1; 0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
units_str_comb={'No signal','Cloud liquid','Drizzle','Rain',...
    'SLW','Ice crystals','Snow','Wet snow/rimed ice','Aerosols'};

f4=figure('DefaultAxesFontSize',12,'Position',[400 300 1300 500]);

s1=subplot(1,1,1);
fig1=surf(time,asl,PID+1,'edgecolor','none');
view(2);
ylim(ylimits);
%xlim([data.time(1),data.time(end)]);
xlim([datetime(2018,1,29,1,50,0),time(end)]);
caxis([.5 9.5]);
colormap(s1,cscale_comb);
cb=colorbar;
cb.Ticks=1:9;
cb.TickLabels=units_str_comb;
ylabel('Altitude (km)');
title(['Particle ID Combined']);