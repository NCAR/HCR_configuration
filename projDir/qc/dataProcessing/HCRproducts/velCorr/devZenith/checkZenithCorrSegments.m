% Plot up and down pointing segments for the zenith pointing corrections

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2
qcVersion='v3.1';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

flight=7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/v3.0_full/10hz/';

[~,indirCorr]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'velCorrZenithPlots/segments/'];

infileZ=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velCorrZenithSegments_',project,'_RF',num2str(flight),'.txt'];
caseListZ = table2array(readtable(infileZ));

infileN=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velCorrNadirSegments_',project,'_RF',num2str(flight),'.txt'];
caseListN = table2array(readtable(infileN));

%% Load vel
disp('Loading velCorr ...');

fileIn1=dir([indirCorr,whichModel,'.velCorr.*.Flight',num2str(flight),'.mat']);
velCorrIn=load([indirCorr,fileIn1.name]);
velCorrGet=velCorrIn.velCorr;

fileIn2=dir([indirCorr,whichModel,'.time.*.Flight',num2str(flight),'.mat']);
timeIn=load([indirCorr,fileIn2.name]);
timeGet=timeIn.timeHCR;

zSegsVEL=[];
zSegsVELcorr=[];
zAsl=[];
nSegsVEL=[];
nAsl=[];

%% Get zenith pointing
for aa=1:size(caseListZ,1)
        
    startTime=datetime(caseListZ(aa,1:6));
    endTime=datetime(caseListZ(aa,7:12));
       
    % Get vel data
    data=[];
    data.VEL_MASKED=[];
   
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    data.VEL_MASKED=-data.VEL_MASKED;

    % Calculate mean vel over segment

    velCorrInds=find(timeGet>=startTime & timeGet<=endTime);
    velCorrThis=velCorrGet(:,velCorrInds);
    velCorrThis(isnan(data.VEL_MASKED))=nan;
    velCorrThis=-velCorrThis;
    dataSum=sum(~isnan(data.VEL_MASKED),2);

    meanVel=mean(data.VEL_MASKED,2,'omitnan');
    meanVel(dataSum<100)=nan;
    meanVelCorr=mean(velCorrThis,2,'omitnan');
    meanVelCorr(dataSum<100)=nan;

    zSegsVEL=cat(2,zSegsVEL,meanVel);
    zSegsVELcorr=cat(2,zSegsVELcorr,meanVelCorr);
    zAsl=cat(2,zAsl,mean(data.asl,2,'omitnan'));
end

%% Get nadir pointing
for aa=1:size(caseListN,1)
        
    startTime=datetime(caseListN(aa,1:6));
    endTime=datetime(caseListN(aa,7:12));
       
    % Get vel data
    data=[];
    data.VEL_MASKED=[];
   
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    % Calculate mean vel over segment
    dataSum=sum(~isnan(data.VEL_MASKED),2);

    meanVel=mean(data.VEL_MASKED,2,'omitnan');
    meanVel(dataSum<100)=nan;
    
    nSegsVEL=cat(2,nSegsVEL,meanVel);   
    nAsl=cat(2,nAsl,mean(data.asl,2,'omitnan'));
end

%% Plot

disp('Plotting ...');
zcols=winter(size(zSegsVEL,2));
ncols=autumn(size(nSegsVEL,2));

close all

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,700],'renderer','painters','visible','on');

% Uncorrected
subplot(1,2,1)
hold on
plot([0,0],[0,10000],'-k','LineWidth',1);
for ii=1:size(zSegsVEL,2)
    plot(zSegsVEL(:,ii),zAsl(:,ii),'-','LineWidth',2,'Color',zcols(ii,:));
end

for ii=1:size(nSegsVEL,2)
    plot(nSegsVEL(:,ii),nAsl(:,ii),'-','LineWidth',2,'Color',ncols(ii,:));
end

xlim([-2 2]);
ylim([0 3000]);
grid on
xlabel('Velocity (m s^{-1}')
ylabel('Altitude (m)')
title([{['RF',num2str(flight)]};{'Blue: uncorrected zenith. Red: nadir.'}])

% Corrected
subplot(1,2,2)
hold on
plot([0,0],[0,10000],'-k','LineWidth',1);
for ii=1:size(zSegsVELcorr,2)
    plot(zSegsVELcorr(:,ii),zAsl(:,ii),'-','LineWidth',2,'Color',zcols(ii,:));
end

for ii=1:size(nSegsVEL,2)
    plot(nSegsVEL(:,ii),nAsl(:,ii),'-','LineWidth',2,'Color',ncols(ii,:));
end

xlim([-2 2]);
ylim([0 3000]);
grid on
xlabel('Velocity (m s^{-1}')
ylabel('Altitude (m)')
title([{['RF',num2str(flight)]};{'Blue: corrected zenith. Red: nadir.'}])

set(gcf,'PaperPositionMode','auto')
print([figdir,'segments_RF',num2str(flight),'.png'],'-dpng','-r0');

