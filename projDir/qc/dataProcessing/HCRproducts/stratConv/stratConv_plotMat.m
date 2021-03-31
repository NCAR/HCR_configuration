% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/cloudPuzzle/',project,'/10hz/'];

outdir=['/run/media/romatsch/RSF0006/rsf/stratConvHCR/',project,'Mat/'];

%figdir=['/h/eol/romatsch/hcrCalib/clouds/stratConv/'];
figdir=['/home/romatsch/plots/HCR/stratConv/',project,'/reflStratConv/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=4:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
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
    
    fileIn1=dir([outdir,'era5.stratConv.*.Flight',num2str(aa),'.mat']);
    sc2D=load([outdir,fileIn1.name]);    
    stratConv=sc2D.stratConv;
    
    fileIn2=dir([outdir,'era5.stratConv1D.*.Flight',num2str(aa),'.mat']);
    sc1DIn=load([outdir,fileIn2.name]);    
    stratConv1D=sc1DIn.stratConv1D;
    
    fileIn3=dir([outdir,'era5.time.*.Flight',num2str(aa),'.mat']);
    scTime=load([outdir,fileIn3.name]);
    timeSC=scTime.timeHCR;
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
    
    colMapSC=[1,0,0;
        1,0,1;
        1,0.5,0;
        0.5,0,1;
        0,0,1;
        0,0.7,1;
        0,0,0];
    
    data.DBZ(data.FLAG>1)=nan;
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        
        if sum(sum(~isnan(dbzPlot)))~=0
            aslPlot=data.asl(:,timeInds);
            
            timeIndsSC=find(timeSC>=startPlot & timeSC<=endPlot);
            sc1D=stratConv1D(timeIndsSC);
            
            stratConvPlot=stratConv(:,timeIndsSC);
            stratConvPlot(stratConvPlot==20)=14;
            stratConvPlot(stratConvPlot==21)=15;
            stratConvPlot(stratConvPlot==30)=16;
            
            
            f1 = figure('Position',[200 500 1500 600],'DefaultAxesFontSize',12);
            
            s1=subplot(2,1,1);
            
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
            
            s5=subplot(30,1,30);
            
            hold on
            scat1=scatter(timePlot,ones(size(timePlot)),10,sc1D,'filled');
            set(gca,'clim',[0,1]);
            set(gca,'YTickLabel',[]);
            %s5.Colormap=;
            xlim([timePlot(1),timePlot(end)]);
            s5pos=s5.Position;
            s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];
            
            s4=subplot(2,1,2);
            
            hold on
            surf(timePlot,aslPlot./1000,stratConvPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([0 10]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            s4.Colormap=colMapSC;
            caxis([9.5 16.5]);
            cb=colorbar;
            cb.Ticks=10:16;
            cb.TickLabels={'Isolated tethered conv.','Embedded tethered conv.','Isolated elevated conv.',...
                'Embedded elevated conv.','Stratiform','Strat. with emb. conv.',...
                'Small'};
            set(gca,'XTickLabel',[]);
            grid on
            title('Stratiform/convective partitioning')
            s4pos=s4.Position;
            s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_Flight',num2str(aa),'_stratConv_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        startPlot=endPlot;
    end
end