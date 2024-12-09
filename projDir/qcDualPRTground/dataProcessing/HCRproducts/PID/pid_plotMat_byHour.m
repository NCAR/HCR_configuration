% Plot HCR pid from mat file in hourly plots
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='meow'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc0, qc1, or qc2
qcVersion='v1.0';
dataFreq='10hz_combined';
whichModel='hrrr';

ylimUpper=15;

showPlot='off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,dataFreq);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,dataFreq);

figdir=[indir(1:end-14),'pidPlots/wholeIOPs/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];

units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
    'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

for aa=1:size(caseList,1)
    disp(['IOP ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ=[];
            
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
    
    %% Load pid
    disp('Loading PID.');
    
    fileIn1=dir([modeldir,whichModel,'.pid.*.IOP',num2str(aa),'.mat']);
    pid2D=load([modeldir,fileIn1.name]);    
    pid=pid2D.pid;
        
    fileIn2=dir([modeldir,whichModel,'.time.*.Iop',num2str(aa),'.mat']);
    pidTime=load([modeldir,fileIn2.name]);
    timePID=pidTime.timeHCR;
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
        
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(90);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        timeInds=timeInds(1:3:length(timeInds));
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        
        if sum(sum(~isnan(dbzPlot)))~=0 & size(dbzPlot,2)>10
            aslPlot=data.asl(:,timeInds);
            
            timeIndsPID=find(timePID>=startPlot & timePID<=endPlot);
            timeIndsPID=timeIndsPID(1:3:length(timeIndsPID));
                        
            pidPlot=pid(:,timeIndsPID);
                        
            
            f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'visible','off');
            
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
            box on
            title('Reflectivity (dBZ)')
            s1pos=s1.Position;
            s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
            
            s2=subplot(2,1,2);
            
            hold on
            surf(timePlot,aslPlot./1000,pidPlot,'edgecolor','none');
            view(2);
            colormap(s2,cscale_hcr);
            cb=colorbar;
            cb.Ticks=1:11;
            cb.TickLabels=units_str_hcr;
            ylabel('Altitude (km)');
            title(['HCR particle ID']);
            caxis([.5 11.5]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            
            grid on
            box on
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_IOP',num2str(aa),'_pid_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        startPlot=endPlot;
    end
end