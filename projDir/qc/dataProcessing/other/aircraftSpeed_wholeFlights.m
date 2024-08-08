% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

showPlot='on';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'aircraftSpeed/wholeFlights/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

% Loop through cases

for aa=19:size(caseList,1)

    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    %% Load data

    disp('Loading data ...');

    data=[];

    data.altitude=[];
    data.eastward_velocity=[];
    data.northward_velocity=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    velAircraft=sqrt(data.eastward_velocity.^2+data.northward_velocity.^2);

    %% Plot in hourly increments

    close all

    disp('Plotting ...');
    ylimitsSpeed=[0,280];
    ylimitsAlt=[0,15];

    f1 = figure('Position',[200 500 1300 800],'DefaultAxesFontSize',12,'visible',showPlot);
    t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

    ax=nexttile(1);

    plot(data.time,velAircraft,'-b','LineWidth',2);
    xlim([data.time(1),data.time(end)]);
    ylim(ylimitsSpeed);
    ylabel('Aircraft speed (m s^{-1})')

    yticks(0:10:300);

    grid on
    box on

    title(['Flight ',num2str(aa),'. Aircraft speed']);

    ax=nexttile(2);

    plot(data.time,data.altitude./1000,'-b','LineWidth',2);
    xlim([data.time(1),data.time(end)]);
    ylim(ylimitsAlt);
    ylabel('Aircraft altitude (km)')

    grid on
    box on

    title(['Aircraft altitude']);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_aircraftSpeed_Flight',num2str(aa),'.png'],'-dpng','-r0');

end