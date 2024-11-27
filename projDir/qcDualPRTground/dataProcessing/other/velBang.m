% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc0';
freqData='10hz_combined';
freqData2='100hz_long';
qcVersion='';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);
indir2=HCRdir(project,quality,qcVersion,freqData2);

figdir='/scr/virga1/rsfdata/projects/meow/hcr/cfradial/moments/velBang/';


%% Run processing

% Go through flights
for ii=12:size(caseList,1)

    disp(['Case ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
       
    %% Load cfradial data 10 Hz
    data=[];
    data.VEL_long=[];
    data.VEL_RAW_long=[];
    data.eastward_velocity=[];
    data.northward_velocity=[];
    data.vertical_velocity=[];
    data.azimuth=[];

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    %% Load cfradial data 50 Hz
    data2=[];
    data2.VEL=[];
    data2.VEL_RAW=[];
    % Make list of files within the specified time frame
    fileList2=makeFileList(indir2,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data2=read_HCR(fileList2,data2,startTime,endTime);
  
    %% Calc corrections
    xcorr=sind(data.azimuth).*cosd(data.elevation).*data.eastward_velocity;
    ycorr=cosd(data.azimuth).*cosd(data.elevation).*data.northward_velocity;
    zcorr=sind(data.elevation).*data.vertical_velocity;

    %% Plot DBZ

    close all

    f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12);

    t = tiledlayout(5,1,'TileSpacing','tight','Padding','compact');

    s1=nexttile(1);

    hold on
    scatter(data.time,data.VEL_RAW_long(15,:));
    scatter(data2.time,data2.VEL_RAW(15,:));
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('VEL RAW in radar bang')
    ylim([-2,2]);
    xlim([data.time(1),data.time(end)]);    

    legend('10 Hz long','50 Hz long')

    s2=nexttile(2);

    hold on
    scatter(data.time,data.VEL_long(15,:));
    scatter(data2.time,data2.VEL(15,:));
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('VEL in radar bang')
    ylim([-2,2]);
    xlim([data.time(1),data.time(end)]);    

    legend('10 Hz long','50 Hz long')

    s3=nexttile(3);

    hold on
    scatter(data.time,data.elevation);
    ylim([89.5,90.5]);
    ylabel('Elevation (deg)');

    yyaxis right
    scatter(data.time,data.azimuth);
    view(2);
    ylabel('Azimuth (deg)');
    grid on
    box on
    title('Angles')
    ylim([-5,365]);
    xlim([data.time(1),data.time(end)]);    

    legend('Elevation','Azimuth')

    s4=nexttile(4);

    hold on
    scatter(data.time,data.eastward_velocity);
    scatter(data.time,data.northward_velocity);
    scatter(data.time,data.vertical_velocity);
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('Ground speed')
    ylim([-1,1]);
    xlim([data.time(1),data.time(end)]);    

    legend('East vel','North vel','Vert vel')

    s5=nexttile(5);

    hold on
    scatter(data.time,xcorr);
    scatter(data.time,ycorr);
    scatter(data.time,zcorr);
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('Velocity corrections')
    ylim([-1,1]);
    xlim([data.time(1),data.time(end)]);    

    legend('Xcorr','Ycorr','Zcorr')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velBang_IOP',num2str(ii)],'-dpng','-r0')

end
