% Remove stripes in HX

clear all;
close all;

startTime=datetime(2021,6,25,19,30,0);
endTime=datetime(2021,6,25,20,10,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset
quality='field'; %field, qc1, or qc2
freqData='100hz'; % 10hz, 100hz, or 2hz
qcVersion='';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.1_full/hxPlots/'];

dataDir=HCRdir(project,quality,qcVersion,freqData);


% Make list of files within the specified time frame
fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

PsVoltage=[];
time=[];

for ii=1:length(fileList)
    infile=fileList{ii};

    % Load data
    PsVoltageFile=varFromCfRadialString(infile,'PsVoltage');
    PsVoltage=[PsVoltage,PsVoltageFile];

    timeIn=datetime(str2num(infile(end-48:end-45)),str2num(infile(end-44:end-43)),str2num(infile(end-42:end-41)),...
        str2num(infile(end-39:end-38)),str2num(infile(end-37:end-36)),str2num(infile(end-35:end-34)));

    time=[time,timeIn];

end
%% Plot
ylimits=[-0.5 10];

close all

figure('DefaultAxesFontSize',11,'position',[1,100,1200,600],'renderer','painters');

plot(time,PsVoltage,'-k','LineWidth',1);
grid on
ylabel('PsVoltage')

formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_PsVoltage'],...
    '-dpng','-r0');
