clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/virga1/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.2_full/specParams/specPaperFigs/';

load([figdir,'widthNoreaster.mat']);
load([figdir,'forReviewers/specData/noreaster/noreaster_spec_20150202_125100_to_20150202_125700.mat']);

%% Match data
startTime=datetime(2015,2,2,12,55,0);
endTime=momentsTime.time(end);

timeInds=find(momentsTime.time>=startTime & momentsTime.time<=endTime);
specInds=find(momentsSpecParams.time>=startTime & momentsSpecParams.time<=endTime);

dbz=momentsTime.dbz(:,timeInds);
vel=momentsTime.vel(:,timeInds);
snr=momentsTime.snr(:,timeInds);
asl=momentsTime.asl(:,timeInds);

widthCut=momentsTime.widthCorr(:,timeInds);
widthCorr=momentsSpecParams.width(:,specInds);

cutInds=find(widthCut==0.1 & asl<7000);
nonCutInds=find(widthCut>0.1 & asl<7000);

dbzC=dbz(cutInds);
velC=vel(cutInds);
snrC=snr(cutInds);

dbzNC=dbz(nonCutInds);
velNC=vel(nonCutInds);
snrNC=snr(nonCutInds);

edges.dbz=-15:2:15;
edges.vel=-5:0.4:2;
edges.snr=0:5:80;

dbzHc=histcounts(dbzC,edges.dbz);
dbzHnc=histcounts(dbzNC,edges.dbz);
velHc=histcounts(velC,edges.vel);
velHnc=histcounts(velNC,edges.vel);
snrHc=histcounts(snrC,edges.snr);
snrHnc=histcounts(snrNC,edges.snr);
%% Figure

close all

% Figure
f1 = figure('Position',[200 500 900 300],'DefaultAxesFontSize',12);

t = tiledlayout(1,3,'TileSpacing','tight','Padding','compact');

s1=nexttile(1);
bar(edges.dbz(1:end-1)+(edges.dbz(2)-edges.dbz(1))/2,cat(1,dbzHc./sum(dbzHc).*100,dbzHnc./sum(dbzHnc).*100),1,GroupWidth=1);
xlim([-13,15]);
ylim([0 22]);
xticks(-15:3:15);
xtickangle(0);

xlabel('dBZ')
ylabel('Percent (%)')

title('Reflectivity')

s2=nexttile(2);
bar(edges.snr(1:end-1)+(edges.snr(2)-edges.snr(1))/2,cat(1,snrHc./sum(snrHc).*100,snrHnc./sum(snrHnc).*100),1,GroupWidth=1);
xlim([0,75]);
ylim([0 20]);
xticks(0:10:80);
xtickangle(0);

xlabel('dBZ')
%ylabel('Percent (%)')

leg=legend('Width=0.1 m/s','Width>0.1 m/s','Location','northeast');
leg.ItemTokenSize=[10,15];

title('Signal to noise ratio')

s3=nexttile(3);
bar(edges.vel(1:end-1)+(edges.vel(2)-edges.vel(1))/2,cat(1,velHc./sum(velHc).*100,velHnc./sum(velHnc).*100),1,GroupWidth=1);
xlim([-5,1.5]);
ylim([0 38]);
xticks(-15:1:15);
xtickangle(0);

xlabel('m s^{-1}')
%ylabel('Percent (%)')

title('Velocity')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'cutOff.png'],'-dpng','-r0');
