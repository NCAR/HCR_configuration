% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc0, qc1, or qc2
qcVersion='v3.2';
dataFreq='10hz';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

showPlot='off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,dataFreq);

[modelNC modeldir]=modelDir(project,whichModel,quality,qcVersion,dataFreq);

figdir=[dataDir(1:end-5),'flagPlots/wholeFlights/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=3:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ=[];
    data.TOPO=[];
            
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
    
    %% Load model
    
    fileIn1=dir([modeldir,whichModel,'.flagfield.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    flag=modIn.flagField;
       
    fileIn1=dir([modeldir,whichModel,'.antstat.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    antstat=modIn.antennaStat;
        
    fileIn1=dir([modeldir,whichModel,'.time.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    timeMod=modIn.timeHCR;
    
    %% Plot in hourly increments
    
    disp('Plotting ...');

    ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
        'Bang (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
        'NS cal (10)';'Missing (11)'};

    colMask=[0.4,0.8,1;
        0,0,0;
        0.5,0,0.5;
        0,1,0;
        0.2,0.6,0.2;
        1,0,0;
        0,0,0.6;
        0.7065,0.4415,0.2812;
        0.5,0.5,0.5;
        0.9290,0.8940,0.1250;
        1,0.6,0];
    
    startPlot=startTime;
    emptyInds=0;
    
    while startPlot<endTime

        close all

        endPlot=startPlot+minutes(60);
        timeIndsAll=find(data.time>=startPlot & data.time<=endPlot);

        if isempty(timeIndsAll)
            startPlot=endPlot;
            continue
        end

        getInd=1:20:length(timeIndsAll);
        timeInds=timeIndsAll(getInd);
                
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        topoPlot=data.TOPO(timeInds);
        
        aslPlot=data.asl(:,timeInds);
        
        timeIndsModAll=find(timeMod>=startPlot & timeMod<=endPlot);
        getIndMod=1:20:length(timeIndsModAll);
        timeIndsMod=timeIndsModAll(getIndMod);
        antstatPlot=antstat(timeIndsMod);
                
        flagPlot=flag(:,timeIndsMod);
        
        dbzMaskedPlot=dbzPlot;
        dbzMaskedPlot(flagPlot>1)=nan;
                        
        f1 = figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible',showPlot);
        
        s1=subplot(4,1,1);
        
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
        box on
        title('Reflectivity (dBZ)')
        s1pos=s1.Position;
                      
        s2=subplot(4,1,2);
        
        hold on
        surf(timePlot,aslPlot./1000,flagPlot,'edgecolor','none');
        view(2);
        caxis([1 11]);
        ylabel('Altitude (km)');
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        grid on
        box on
        title('FLAG')
        plot(timePlot,topoPlot./1000,'-k','linewidth',2);
        
        s2.SortMethod = 'childorder';
        
        s2.Colormap=colMask;
        hcb=colorbar;
        hcb.Ticks=[1.5:0.9:11.5];
        hcb.TickLabels=ytickLabels;
               
        s3=subplot(4,1,3);
        
        hold on
        surf(timePlot,aslPlot./1000,dbzMaskedPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-35 25]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        box on
        title('Masked reflectivity (dBZ)')

        s4=subplot(4,1,4);

        hold on
        plot(timePlot,antstatPlot,'-b','linewidth',2);
        xlim([timePlot(1),timePlot(end)]);
        grid on
        box on

        ylabel('Ant. Stat.');
        yticks(1:6)
        yticklabels({'Down (1)','Up (2)','Pointing (3)','Scanning (4)','Transition (5)','Failure(6)'})
        ylim([0 7])

        title('Antenna status')

        s4pos=s4.Position;
        s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
                
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_flag_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        startPlot=endPlot;
    end
end