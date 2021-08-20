% Plot HCR convStrat from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
% dataFreq='10hz';
% qcVersion='v2.1';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

if strcmp(project,'otrec')
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/convStrat/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/convStrat/10hz/';
end

figdir=[indir(1:end-5),'convStratPlots/wholeFlights/'];

modeldir=[indir(1:end-47),'mat/convStrat/10hz/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

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
    
    data.DBZ = [];
    data.FLAG=[];
        
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
       
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    %% Load stratconv
    
    disp('Loading conv/strat data.');
    
    fileIn1=dir([modeldir,'era5.convStrat.*.Flight',num2str(aa),'.mat']);
    sc2D=load([modeldir,fileIn1.name]);    
    stratConv=sc2D.convStrat;
    
    fileIn2=dir([modeldir,'era5.convectivity.*.Flight',num2str(aa),'.mat']);
    convIn=load([modeldir,fileIn2.name]);    
    convectivity=convIn.convectivity;
    
    fileIn3=dir([modeldir,'era5.convStrat1D.*.Flight',num2str(aa),'.mat']);
    sc1DIn=load([modeldir,fileIn3.name]);    
    stratConv1D=sc1DIn.convStrat1D;
    
    fileIn4=dir([modeldir,'era5.time.*.Flight',num2str(aa),'.mat']);
    scTime=load([modeldir,fileIn4.name]);
    timeSC=scTime.timeHCR;
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
    
    data.DBZ(data.FLAG>1)=nan;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        
        if sum(sum(~isnan(dbzPlot)))==0
            startPlot=endPlot;
            continue
        end
        
        aslPlot=data.asl(:,timeInds);
        
        timeIndsSC=find(timeSC>=startPlot & timeSC<=endPlot);
        sc1D=stratConv1D(timeIndsSC);
        sc1D(sc1D==14)=1;
        sc1D(sc1D==16)=2;
        sc1D(sc1D==18)=3;
        sc1D(sc1D==25)=4;
        sc1D(sc1D==30)=5;
        sc1D(sc1D==32)=6;
        sc1D(sc1D==34)=7;
        sc1D(sc1D==36)=8;
        sc1D(sc1D==38)=9;
        
        convPlot=convectivity(:,timeIndsSC);
        
        time1D=timePlot(~isnan(sc1D));
        sc1D=sc1D(~isnan(sc1D));
        col1D=colmapSC(sc1D,:);
        
        stratConvPlot=stratConv(:,timeIndsSC);
        stratConvPlot(stratConvPlot==14)=1;
        stratConvPlot(stratConvPlot==16)=2;
        stratConvPlot(stratConvPlot==18)=3;
        stratConvPlot(stratConvPlot==25)=4;
        stratConvPlot(stratConvPlot==30)=5;
        stratConvPlot(stratConvPlot==32)=6;
        stratConvPlot(stratConvPlot==34)=7;
        stratConvPlot(stratConvPlot==36)=8;
        stratConvPlot(stratConvPlot==38)=9;
        
        if sum(sum(~isnan(stratConvPlot)))==0
            startPlot=endPlot;
            continue
        end
        
        f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12,'visible','off');
        
        s1=subplot(3,1,1);
        
        colormap jet
        
        hold on
        surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-35 25]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        title('Reflectivity (dBZ)')
        s1pos=s1.Position;
        
        s2=subplot(3,1,2);
        
        colormap jet
        
        hold on
        surf(timePlot,aslPlot./1000,convPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([0 1]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        title('Convectivity');
        s2pos=s2.Position;
        s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
        
        s4=subplot(30,1,30);
        
        hold on
        scat1=scatter(time1D,ones(size(time1D)),10,col1D,'filled');
        %set(gca,'clim',[0,1]);
        set(gca,'YTickLabel',[]);
        s4.Colormap=colmapSC;
        xlim([timePlot(1),timePlot(end)]);
        s4pos=s4.Position;
        s4.Position=[s4pos(1),s4pos(2)-0.023,s1pos(3),s4pos(4)];
        
        s3=subplot(3,1,3);
        
        hold on
        surf(timePlot,aslPlot./1000,stratConvPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([0 10]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        s3.Colormap=colmapSC;
        caxis([0.5 9.5]);
        cb=colorbar;
        cb.Ticks=1:9;
        cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
            'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
        set(gca,'XTickLabel',[]);
        grid on
        title('Stratiform/convective partitioning')
        s3pos=s3.Position;
        s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_convStrat_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        startPlot=endPlot;
    end
end