% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc0';
freqData='100hz_long';
qcVersion='';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/noise_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir='/scr/virga1/rsfdata/projects/meow/hcr/cfradial/moments/noise/';


%% Run processing

% Go through flights
for ii=2:size(caseList,1)

    disp(['Case ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBMVC=[];
    data.DBMHX=[];
   
    %% Load cfradial data
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    asl=HCRrange2asl(data.range,data.elevation,data.altitude);

     %% Load xml data

    % BodyCurrent=[];
    % CollectorCurrent=[];
    % CathodeVoltage=[];
    % Locked15_5GHzPLO=[];
    % timeTemp=[];
    % 
    % for jj=1:length(fileList)
    %     tempFile=fileList{jj};
    %     BodyCurrent=[BodyCurrent varFromCfRadialString(tempFile,'BodyCurrent')];
    %     CollectorCurrent=[CollectorCurrent varFromCfRadialString(tempFile,'CollectorCurrent')];
    %     CathodeVoltage=[CathodeVoltage varFromCfRadialString(tempFile,'CathodeVoltage')];
    %     Locked15_5GHzPLO=[Locked15_5GHzPLO varFromCfRadialString(tempFile,'Locked15_5GHzPLO')];
    % 
    %     startTimeIn=ncread(tempFile,'time_coverage_start')';
    %     startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
    %         str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    %     timeTemp=[timeTemp startTimeFile];
    % end

    tempFile='/scr/virga1/rsfdata/projects/meow/hcr/qc0/txt/MEOW.temperatures.txt';
    tempnames={'count','year','month','day','hour','min','sec','unix_time',...
        'unix_day','XmitterTemp','PloTemp','EikTemp','VLnaTemp','HLnaTemp',...
        'PolarizationSwitchTemp','RfDetectorTemp','NoiseSourceTemp','Ps28VTemp',...
        'RdsInDuctTemp','RotationMotorTemp','TiltMotorTemp','CmigitsTemp',...
        'TailconeTemp','PentekFpgaTemp','PentekBoardTemp', ...
        'BodyCurrent','CathodeVoltage','CollectorCurrent','Locked15_5GHzPLO'};
    indata=readtable(tempFile);
    indata.Properties.VariableNames=tempnames;

    timeTemp=datetime(indata.year,indata.month,indata.day,indata.hour,indata.min,indata.sec);

    %% Calculate noise power of last gates

    meanVC=mean(data.DBMVC(749:end,:),1);
    meanHX=mean(data.DBMHX(749:end,:),1);

    %% Plot DBZ

    close all

    ylims=[0,12];
    climsDBZ=[-50,25];
    pix=100;
    colJet=cat(1,[0,0,0],jet);
    
    f1 = figure('Position',[200 500 1500 1500],'DefaultAxesFontSize',12);

    t = tiledlayout(3,1,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    hold on
    surf(data.time(1:pix:length(data.time)),data.range(:,1:pix:length(data.time))./1000,data.DBMVC(:,1:pix:length(data.time)),'edgecolor','none');
    view(2);
    ylabel('Range (km)');
    clim([-110,-75]);
    s1.Colormap=jet;
    colorbar
    grid on
    box on
    title('DBMVC long (dB)')
    ylim(ylims);
    xlim([data.time(1),data.time(end)]);    

    s2=nexttile(2);

    hold on
    plot(data.time,meanVC,'-b','LineWidth',1.5);
    plot(data.time,meanHX,'-r','LineWidth',1.5);
    ylabel('Power (dB)');

    legend('DBMVC','DBMHX')

    grid on
    box on
    title('Mean power of last 10 gates (dB)')
    %ylim(ylims);
    xlim([data.time(1),data.time(end)]);

    s3=nexttile(3);

    hold on
    plot(timeTemp,indata.BodyCurrent,'-','LineWidth',1.5);
    plot(timeTemp,indata.CollectorCurrent,'-','LineWidth',1.5);
    plot(timeTemp,indata.Locked15_5GHzPLO,'-','LineWidth',1.5);

    yyaxis right
    plot(timeTemp,indata.CathodeVoltage,'-','LineWidth',1.5);

    legend({'BodyCurrent','CollectorCurrent','Locked15_5GhzPLO','CathodeVoltage'}, ...
        'Location','southoutside','Orientation','horizontal','Interpreter','none');

    grid on
    box on
    title('Status')
    %ylim(ylims);
    xlim([data.time(1),data.time(end)]);

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_noise_',datestr(data.time(1),'yyyymmdd_HHMM')],'-dpng','-r0')

end
