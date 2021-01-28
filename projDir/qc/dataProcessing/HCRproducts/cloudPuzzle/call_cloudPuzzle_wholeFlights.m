% Analyze HCR clouds

clear all;
close all;

plotTest=1;

project='cset'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

if strcmp(project,'otrec')
    ylimits=[-0.2 15];
else
    ylimits=[-0.2 10];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudPuzzle/'];
figdir=['/home/romatsch/plots/HCR/cloudPuzzle/',project,'/plotAll/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

outdir=['/run/media/romatsch/RSF0006/rsf/cloudPuzzle/',project,'Mat/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
        
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
        
    %% Load data
        
    data=[];
    
    data.DBZ=[];
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
     %% Cloud puzzle
    disp('Making cloud puzzle');
    
    cloudPuzzle=f_cloudPuzzle_radial(data);
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
    
    while startPlot<endTime
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        aslPlot=data.asl(:,timeInds);
        puzzlePlot=cloudPuzzle(:,timeInds);
        
        uClouds=unique(puzzlePlot(~isnan(puzzlePlot)));
        if isempty(uClouds)
            startPlot=endPlot;
            continue
        end
        
        if uClouds(1)~=0
            uClouds=[0;uClouds];
        end
        cloudCount=length(uClouds);
        
        if cloudCount>0 & length(uClouds)~=1
            minCloud=uClouds(2);
            maxCloud=uClouds(end);
            
            puzzlePlot(puzzlePlot==0)=minCloud-1;
            
            close all
            
            fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
            
            %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ax1=subplot(2,1,1);
            hold on;
            
            sub1=surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
            view(2);
            sub1=colMapDBZ(sub1);
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([timePlot(1),timePlot(end)]);
            title('Reflectivity')
            grid on
            
            ax2=subplot(2,1,2);
            
            colMapIn=jet(maxCloud-minCloud+1);
            % Make order random
            indsCol=randperm(size(colMapIn,1));
            colMapInds=cat(2,indsCol',colMapIn);
            colMapInds=sortrows(colMapInds);
            colMap=cat(1,[0 0 0],colMapInds(:,2:end));
            
            hold on;
            sub2=surf(timePlot,aslPlot./1000,puzzlePlot,'edgecolor','none');
            view(2);
            ax2.Colormap=colMap;
            caxis([minCloud-1.5 maxCloud+0.5])
            ylim(ylimits);
            ylabel('Altitude (km)');
            xlim([timePlot(1),timePlot(end)]);
            title('Cloud Puzzle')
            grid on
            cb=colorbar;
            cb.TickLabels={''};
            
            formatOut = 'yyyymmdd_HHMM';
            set(gcf,'PaperPositionMode','auto')
            print([figdir,project,'_Flight',num2str(aa),'_',datestr(timePlot(1),formatOut),'_to_',datestr(timePlot(end),formatOut),'_cloudPuzzle.png'],'-dpng','-r0');
        end
        startPlot=endPlot;
    end
    
    %% Saving data
    disp('Saving meltLayer field ...')
    
    save([outdir,whichModel,'.cloudPuzzle.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'cloudPuzzle');
    
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
end
