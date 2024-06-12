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

figdir=[indir(1:end-14),'flightPlots/'];


%% Run processing

% Go through flights
for ii=7:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_long=[];
   
    %% Load data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    asl=HCRrange2asl(data.range,data.elevation,data.altitude);

    %% Plot DBZ

    close all

    ylims=[0,12];
    climsDBZ=[-50,25];
    pix=100;
    colJet=cat(1,[0,0,0],jet);
    
    f1 = figure('Position',[200 500 1500 700],'DefaultAxesFontSize',12);

    t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    dbzTemp=data.DBZ_long;
    dbzTemp(isnan(dbzTemp))=-999;

    hold on
    surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,dbzTemp(:,1:pix:length(data.time)),'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    clim(climsDBZ);
    s1.Colormap=colJet;
    colorbar
    grid on
    box on
    title('DBZ long (dBZ)')
    ylim(ylims);
    xlim([data.time(1),data.time(end)]);    

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_DBZ_IOP',num2str(ii)],'-dpng','-r0')

end
