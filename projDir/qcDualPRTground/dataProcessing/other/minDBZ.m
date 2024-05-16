% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc0';
freqData='10hz_combined';
qcVersion='';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'sensitivity/'];

edges=-100:0.5:100;
edgesS=-30:0.5:100;

dbzAllS=zeros(1,length(edges)-1);
dbz1kmS=zeros(1,length(edges)-1);
dbz6kmS=zeros(1,length(edges)-1);

snrAllS=zeros(1,length(edgesS)-1);
snr1kmS=zeros(1,length(edgesS)-1);
snr6kmS=zeros(1,length(edgesS)-1);

dbzAllL=zeros(1,length(edges)-1);
dbz1kmL=zeros(1,length(edges)-1);
dbz6kmL=zeros(1,length(edges)-1);

snrAllL=zeros(1,length(edgesS)-1);
snr1kmL=zeros(1,length(edgesS)-1);
snr6kmL=zeros(1,length(edgesS)-1);

%% Run processing

% Go through flights
for ii=2:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_short=[];
    data.DBZ_long=[];
    data.SNRVC_short=[];
    data.SNRVC_long=[];

    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    asl=HCRrange2asl(data.range,data.elevation,data.altitude);

    dbzStemp=data.DBZ_short;
    dbzLtemp=data.DBZ_long;
    snrStemp=data.SNRVC_short;
    snrLtemp=data.SNRVC_long;

    % Censor data
    maskS=~isnan(dbzStemp);
    maskL=~isnan(dbzLtemp);
    maskS(1:30,:)=0;
    maskL(1:30,:)=0;
    maskS(isinf(dbzStemp))=0;
    maskL(isinf(dbzLtemp))=0;

    maskS=bwareaopen(maskS,200);
    maskL=bwareaopen(maskL,200);

    dbzStemp(maskS==0)=nan;
    dbzLtemp(maskL==0)=nan;
    snrStemp(maskS==0)=nan;
    snrLtemp(maskL==0)=nan;

    histDBZs=histcounts(dbzStemp,edges);
    dbzAllS=dbzAllS+histDBZs;
    histDBZl=histcounts(dbzLtemp,edges);
    dbzAllL=dbzAllL+histDBZl;

    histSNRs=histcounts(snrStemp,edgesS);
    snrAllS=snrAllS+histSNRs;
    histSNRl=histcounts(snrLtemp,edgesS);
    snrAllL=snrAllL+histSNRl;

    range1km=find(data.range>=800 & data.range <=1200);
    range6km=find(data.range>=5900 & data.range <=6100);

    histDBZ1kmS=histcounts(dbzStemp(range1km),edges);
    histDBZ6kmS=histcounts(dbzStemp(range6km),edges);
    histDBZ1kmL=histcounts(dbzLtemp(range1km),edges);
    histDBZ6kmL=histcounts(dbzLtemp(range6km),edges);

    histSNR2kmS=histcounts(snrStemp(range1km),edgesS);
    histSNR8kmS=histcounts(snrStemp(range6km),edgesS);
    histSNR2kmL=histcounts(snrLtemp(range1km),edgesS);
    histSNR8kmL=histcounts(snrLtemp(range6km),edgesS);

    dbz1kmS=dbz1kmS+histDBZ1kmS;
    dbz6kmS=dbz6kmS+histDBZ6kmS;
    dbz1kmL=dbz1kmL+histDBZ1kmL;
    dbz6kmL=dbz6kmL+histDBZ6kmL;

    snr1kmS=snr1kmS+histSNR2kmS;
    snr6kmS=snr6kmS+histSNR8kmS;
    snr1kmL=snr1kmL+histSNR2kmL;
    snr6kmL=snr6kmL+histSNR8kmL;

    shortLongFrac=sum(maskS(:))/sum(maskL(:));

end
%% Bar plot
close all
f1 = figure('Position',[200 500 2500 1250],'DefaultAxesFontSize',12);

t = tiledlayout(3,4,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

b=bar(edges(1:end-1)+0.25,cat(1,dbzAllS,dbzAllL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title([project,' DBZ']);
xlim([-65 20]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s2=nexttile(2);
bar(edges(1:end-1)+0.25,cat(1,dbzAllS,dbzAllL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title(['DBZ lower end']);
xlim([-65 -45]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s3=nexttile(3);
bar(edgesS(1:end-1)+0.25,cat(1,snrAllS,snrAllL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC']);
xlim([-20 50]);
grid on
box on
legend({'Short';'Long'},'Location','northeast')

s4=nexttile(4);
bar(edgesS(1:end-1)+0.25,cat(1,snrAllS,snrAllL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC lower end']);
xlim([-20 0]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s5=nexttile(5);
b=bar(edges(1:end-1)+0.25,cat(1,dbz1kmS,dbz1kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title([project,' DBZ 1 km']);
xlim([-55 15]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s6=nexttile(6);
bar(edges(1:end-1)+0.25,cat(1,dbz1kmS,dbz1kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title(['DBZ 1 km lower end']);
xlim([-55 -35]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s7=nexttile(7);
bar(edgesS(1:end-1)+0.25,cat(1,snr1kmS,snr1kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC 1 km']);
xlim([-20 45]);
grid on
box on
legend({'Short';'Long'},'Location','northeast')

s8=nexttile(8);
bar(edgesS(1:end-1)+0.25,cat(1,snr1kmS,snr1kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC 1 km lower end']);
xlim([-20 0]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s9=nexttile(9);
b=bar(edges(1:end-1)+0.25,cat(1,dbz6kmS,dbz6kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title([project,' DBZ 6 km']);
xlim([-40 -10]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s10=nexttile(10);
bar(edges(1:end-1)+0.25,cat(1,dbz6kmS,dbz6kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dBZ');
title(['DBZ 6 km lower end']);
xlim([-40 -20]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s11=nexttile(11);
bar(edgesS(1:end-1)+0.25,cat(1,snr6kmS,snr6kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC 6 km']);
xlim([-20 5]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

s12=nexttile(12);
bar(edgesS(1:end-1)+0.25,cat(1,snr6kmS,snr6kmL),'GroupWidth',1,'BarWidth',1);
xlabel('dB');
title(['SNRVC 6 km lower end']);
xlim([-20 -10]);
grid on
box on
legend({'Short';'Long'},'Location','northwest')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_sensitivity'],'-dpng','-r0')

%% Reflectivity

close all

%colDiff=cat(1,[0,0,0],velCols);
colDiff=cat(1,[0,0,0],jet(21));
colJet=cat(1,[0,0,0],jet);
ylims=[0,8];
climsDBZ=[-50,15];
climsSNR=[-30,50];
climsDiffDBZ=[-5,5];
climsDiffSNR=[-10,10];

pix=100;

f1 = figure('Position',[200 500 2400 1200],'DefaultAxesFontSize',12);

t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

dbzStemp(isnan(dbzStemp))=-999;

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,dbzStemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsDBZ);
s1.Colormap=colJet;
colorbar
grid on
box on
title('DBZ short (dBZ)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

s2=nexttile(2);

snrStemp(isnan(snrStemp))=-999;

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,snrStemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsSNR);
s2.Colormap=colJet;
colorbar
grid on
box on
title('SNR short (dB)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

s3=nexttile(3);

dbzLtemp(isnan(dbzLtemp))=-99;

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,dbzLtemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsDBZ);
s3.Colormap=colJet;
colorbar
grid on
box on
title('DBZ long (dBZ)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

s4=nexttile(4);

snrLtemp(isnan(snrLtemp))=-99;

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,snrLtemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsSNR);
s4.Colormap=colJet;
colorbar
grid on
box on
title('SNR long (dB)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

s5=nexttile(5);

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,dbzStemp(:,1:pix:length(data.time))-dbzLtemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsDiffDBZ);
s5.Colormap=colDiff;
colorbar
grid on
box on
title('DBZ short - long (dB)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

s6=nexttile(6);

hold on
surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,snrStemp(:,1:pix:length(data.time))-snrLtemp(:,1:pix:length(data.time)),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim(climsDiffSNR);
s6.Colormap=colDiff;
colorbar
grid on
box on
title('SNR short - long (dB)')
ylim(ylims);
xlim([data.time(1),data.time(end)]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_reflectivity'],'-dpng','-r0')