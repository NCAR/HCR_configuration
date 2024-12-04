% Analyze vel in bang
clearvars;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc2';
freqData='10hz_combined';
freqData2='50hz_longPulse';
qcVersion='v1.0';
load10hz=0;
load50hz=1;

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);
indir2=HCRdir(project,quality,qcVersion,freqData2);

figdir='/scr/virga1/rsfdata/projects/meow/hcr/cfradial/moments/velBang/';


%% Run processing

% Go through flights
for ii=1:size(caseList,1)

    disp(['Case ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    %% Load cfradial data 10 Hz
    if load10hz
        disp('Loading 10 hz data ...')
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
    end

    %% Load cfradial data 50 Hz
    if load50hz
        disp('Loading 50 hz data ...')
        data2=[];
        data2.VEL=[];
        data2.VEL_RAW=[];
        data2.eastward_velocity=[];
        data2.northward_velocity=[];
        data2.vertical_velocity=[];
        data2.azimuth=[];
        % Make list of files within the specified time frame
        fileList2=makeFileList(indir2,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

        % Load data
        data2=read_HCR(fileList2,data2,startTime,endTime);
    end

    %% Calc corrections
    if load10hz
        xcorr=sind(data.azimuth).*cosd(data.elevation).*data.eastward_velocity;
        ycorr=cosd(data.azimuth).*cosd(data.elevation).*data.northward_velocity;
        zcorr=sind(data.elevation).*data.vertical_velocity;
    else
        xcorr=sind(data2.azimuth).*cosd(data2.elevation).*data2.eastward_velocity;
        ycorr=cosd(data2.azimuth).*cosd(data2.elevation).*data2.northward_velocity;
        zcorr=sind(data2.elevation).*data2.vertical_velocity;
    end

    %% Plot DBZ

    close all

    f1 = figure('Position',[200 500 1000 1100],'DefaultAxesFontSize',12);

    t = tiledlayout(5,1,'TileSpacing','tight','Padding','compact');

    s1=nexttile(1);

    hold on
    if load10hz
        scatter(data.time,data.VEL_RAW_long(15,:));
    end
    if load50hz
        scatter(data2.time,data2.VEL_RAW(15,:));
    end
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('VEL RAW in radar bang')
    ylim([-2,2]);
    if load10hz
        xlim([data.time(1),data.time(end)]);
    else
        xlim([data2.time(1),data2.time(end)]);
    end

    if load10hz & load50hz
        legend('10 Hz long','50 Hz long')
    end

    s2=nexttile(2);

    hold on
    if load10hz
        scatter(data.time,data.VEL_long(15,:));
    end
    if load50hz
        scatter(data2.time,data2.VEL(15,:));
    end
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('VEL in radar bang')
    ylim([-2,2]);
    if load10hz
        xlim([data.time(1),data.time(end)]);
    else
        xlim([data2.time(1),data2.time(end)]);
    end

    if load10hz & load50hz
        legend('10 Hz long','50 Hz long')
    end

    s3=nexttile(3);

    hold on
    if load10hz
        scatter(data.time,data.elevation);
    else
        scatter(data2.time,data2.elevation);
    end
    ylim([89.5,90.5]);
    ylabel('Elevation (deg)');

    yyaxis right
    if load10hz
        scatter(data.time,data.azimuth);
    else
        scatter(data2.time,data2.azimuth);
    end
    view(2);
    ylabel('Azimuth (deg)');
    grid on
    box on
    title('Angles')
    ylim([-5,365]);
    if load10hz
        xlim([data.time(1),data.time(end)]);
    else
        xlim([data2.time(1),data2.time(end)]);
    end

    legend('Elevation','Azimuth')

    s4=nexttile(4);

    hold on
    if load10hz
        scatter(data.time,data.eastward_velocity);
        scatter(data.time,data.northward_velocity);
        scatter(data.time,data.vertical_velocity);
    else
        scatter(data2.time,data2.eastward_velocity);
        scatter(data2.time,data2.northward_velocity);
        scatter(data2.time,data2.vertical_velocity);
    end
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('Ground speed')
    ylim([-1,1]);
    if load10hz
        xlim([data.time(1),data.time(end)]);
    else
        xlim([data2.time(1),data2.time(end)]);
    end

    legend('East vel','North vel','Vert vel')

    s5=nexttile(5);

    hold on
    if load10hz
        scatter(data.time,xcorr);
        scatter(data.time,ycorr);
        scatter(data.time,zcorr);
    else
        scatter(data2.time,xcorr);
        scatter(data2.time,ycorr);
        scatter(data2.time,zcorr);
    end
    view(2);
    ylabel('Velocity (m/s)');
    grid on
    box on
    title('Velocity corrections')
    ylim([-1,1]);
    if load10hz
        xlim([data.time(1),data.time(end)]);
    else
        xlim([data2.time(1),data2.time(end)]);
    end

    legend('Xcorr','Ycorr','Zcorr')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velBang_',quality,'IOP',num2str(ii)],'-dpng','-r0')

end
