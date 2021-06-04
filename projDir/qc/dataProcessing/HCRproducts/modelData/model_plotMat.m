% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset, otrec
quality='qc0'; %field, qc0, qc1, or qc2
qcVersion='v0.1';
dataFreq='10hz';
whichModel='ecmwf';

if strcmp(project,'otrec') | strcmp(project,'spicule')
    ylimUpper=15;
else
    ylimUpper=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,dataFreq);

[modelNC modeldir]=modelDir(project,whichModel,dataFreq);

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc0/cfradial/v0.1/modelPlots/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
            
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
    
    fileIn1=dir([modeldir,whichModel,'.pHCR.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    press=modIn.pHCR;
    
    fileIn1=dir([modeldir,whichModel,'.rhHCR.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    rh=modIn.rhHCR;
    
    fileIn1=dir([modeldir,whichModel,'.tempHCR.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    temperature=modIn.tempHCR;
    
%     fileIn1=dir([modeldir,whichModel,'.sst.*.Flight',num2str(aa),'.mat']);
%     modIn=load([modeldir,fileIn1.name]);    
%     sst=modIn.sstHCR;
    
    fileIn1=dir([modeldir,whichModel,'.topo.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    topo=modIn.topo;
    
    fileIn1=dir([modeldir,whichModel,'.uSurf.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    u=modIn.uSurfHCR;
    
    fileIn1=dir([modeldir,whichModel,'.vSurf.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    v=modIn.vSurfHCR;
    
    fileIn1=dir([modeldir,whichModel,'.time.*.Flight',num2str(aa),'.mat']);
    modIn=load([modeldir,fileIn1.name]);    
    timeMod=modIn.timeHCR;
        
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(60);
        timeIndsAll=find(data.time>=startPlot & data.time<=endPlot);
        getInd=1:20:length(timeIndsAll);
        timeInds=timeIndsAll(getInd);
                
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        
        aslPlot=data.asl(:,timeInds);
        
        timeIndsModAll=find(timeMod>=startPlot & timeMod<=endPlot);
        getIndMod=1:20:length(timeIndsModAll);
        timeIndsMod=timeIndsModAll(getIndMod);
%         sstPlot=sst(timeIndsMod);
        topoPlot=topo(timeIndsMod);
        uPlot=u(timeIndsMod);
        vPlot=v(timeIndsMod);
        
        pressPlot=press(:,timeIndsMod);
        rhPlot=rh(:,timeIndsMod);
        tempPlot=temperature(:,timeIndsMod);
                
        f1 = figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12);
        
        s1=subplot(5,1,1);
        
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
                      
        s2=subplot(5,1,2);
        
        hold on
        surf(timePlot,aslPlot./1000,pressPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([200 1000]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        box on
        title('Pressure (hPa)')
        
        s3=subplot(5,1,3);
        
        hold on
        surf(timePlot,aslPlot./1000,rhPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        s3.Colormap=flipud(jet);
        caxis([0 100]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        box on
        title('Relative humidity (%)')
        
        s4=subplot(5,1,4);
        
        hold on
        surf(timePlot,aslPlot./1000,tempPlot,'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-50 30]);
        ylim([0 ylimUpper]);
        xlim([timePlot(1),timePlot(end)]);
        colorbar
        grid on
        box on
        title('Temperature (C)')
        
        s5=subplot(5,1,5);
        
        hold on
        plot(timePlot,uPlot,'-b','linewidth',2);
        plot(timePlot,vPlot,'-g','linewidth',2);
        ylabel('Wind (m s^{-1})');
        ylim([-15 15]);
        yticks(-15:5:15);
        xlim([timePlot(1),timePlot(end)]);
        grid on
        box on
        
        yyaxis right
        plot(timePlot,topoPlot,'-r','linewidth',2);
        ylabel('Topography (m)');
        ylim([0 3000]);
        yticks(0:500:3000);
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k';
        
        title('Surface wind and topography')
        legend('u','v','topo','orientation','horizontal');
        
        s5pos=s5.Position;
        s5.Position=[s5pos(1) s5pos(2) s1pos(3) s5pos(4)];
                
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_model_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        startPlot=endPlot;
    end
end