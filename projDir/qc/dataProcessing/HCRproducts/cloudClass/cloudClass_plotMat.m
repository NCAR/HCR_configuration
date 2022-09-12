% Plot HCR convStrat from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

showPlot='off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'cloudClass/wholeFlights/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

colmapCC=[0,0,0;
    204,255,204;
    153,204,0;
    0,128,0;
    0,204,255;
    51,102,255;
    0,0,180;
    255,128,128;
    255,0,0;
    150,0,0;
    255,204,0;
    255,153,0;
    255,102,0];

colmapCC=colmapCC./255;

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    0.7,0,0];

for aa=2:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.ECHO_TYPE_2D=[];
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
           
    %% Load cloudClass
    
    disp('Loading cloudClass data.');
    
    fileIn1=dir([modeldir,'era5.cloudClass.*.Flight',num2str(aa),'.mat']);
    cc=load([modeldir,fileIn1.name]);    
    cloudClass=cc.cloudClass;

    fileIn2=dir([modeldir,'era5.time.*.Flight',num2str(aa),'.mat']);
    ccTime=load([modeldir,fileIn2.name]);
    timeCC=ccTime.timeHCR;

    %% Plot in hourly increments

    disp('Plotting ...');

    startPlot=startTime;

    while startPlot<endTime

        close all

        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);

        timePlot=data.time(timeInds);

        echoPlot=data.ECHO_TYPE_2D(:,timeInds);
        echoPlot(echoPlot==14)=1;
        echoPlot(echoPlot==16)=2;
        echoPlot(echoPlot==18)=3;
        echoPlot(echoPlot==25)=4;
        echoPlot(echoPlot==30)=5;
        echoPlot(echoPlot==32)=6;
        echoPlot(echoPlot==34)=7;
        echoPlot(echoPlot==36)=8;
        echoPlot(echoPlot==38)=9;

        if sum(sum(~isnan(echoPlot)))==0
            startPlot=endPlot;
            continue
        end

        aslPlot=data.asl(:,timeInds);

        timeIndsCC=find(timeCC>=startPlot & timeCC<=endPlot);

        classPlot=cloudClass(:,timeIndsCC);
        classPlot(classPlot==11)=4;
        classPlot(classPlot==12)=5;
        classPlot(classPlot==13)=6;
        classPlot(classPlot==21)=7;
        classPlot(classPlot==22)=8;
        classPlot(classPlot==23)=9;
        classPlot(classPlot==31)=10;
        classPlot(classPlot==32)=11;
        classPlot(classPlot==33)=12;

        if sum(sum(~isnan(classPlot)))==0
            startPlot=endPlot;
            continue
        end

        f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible',showPlot);

        s1=subplot(2,1,1);

        hold on
        surf(timePlot,aslPlot./1000,echoPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([0 10]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        s1.Colormap=colmapSC;
        caxis([0.5 9.5]);
        cb=colorbar;
        cb.Ticks=1:9;
        cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
            'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
        grid on
        box on
        title('Echo type')

        s2=subplot(2,1,2);

        hold on
        surf(timePlot,aslPlot./1000,classPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-1 13]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        s2.Colormap=colmapCC;
        caxis([-0.5 12.5]);
        cb=colorbar;
        cb.Ticks=0:12;
        cb.TickLabels={'Not classified','Strat Low','Strat Mid','Strat High',...
            'Strat Precip Shallow','Strat Precip Mid','Strat Precip Deep',...
            'Conv Young Shallow','Conv Young Mid','Conv Yong Deep',...
            'Conv Mature Shallow','Conv Mature Mid','Conv Mature Deep'};
        grid on
        box on
        title('Cloud classification')

        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_cloudClass_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')

        startPlot=endPlot;
    end
end