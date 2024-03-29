% Plot HCR pid from mat file in hourly plots

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
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
end

figdir=[indir(1:end-5),'pidPlots/wholeFlights/temperatures/'];

modeldir=[indir(1:end-30),'mat/pid/10hz/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid','Mixed Phase','Large Frozen','Small Frozen'};
        

for aa=2:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.TEMP = [];
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
    
    %% Load pid
    disp('Loading PID.');
    
    fileIn1=dir([modeldir,whichModel,'.pid.*.Flight',num2str(aa),'.mat']);
    pid2D=load([modeldir,fileIn1.name]);    
    pid=pid2D.pid;
        
    fileIn2=dir([modeldir,whichModel,'.time.*.Flight',num2str(aa),'.mat']);
    pidTime=load([modeldir,fileIn2.name]);
    timePID=pidTime.timeHCR;
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    startPlot=startTime;
           
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        tempPlot=data.TEMP(1,timeInds);
        
        if sum(sum(~isnan(tempPlot)))~=0
            aslPlot=data.asl(:,timeInds);
            
            timeIndsPID=find(timePID>=startPlot & timePID<=endPlot);
                        
            pidPlot=pid(:,timeIndsPID);
                        
            
            f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'visible','off');
            
            s1=subplot(2,1,1);
            
            colormap jet
            
            hold on
            plot(timePlot,tempPlot,'-k','linewidth',2);
            ylabel('Temperature (C)');
            caxis([-35 25]);
            ylim([-70 40]);
            xlim([timePlot(1),timePlot(end)]);
            grid on
            box on
            title('Temperature')
            s1pos=s1.Position;
            
            s2=subplot(2,1,2);
            
            hold on
            surf(timePlot,aslPlot./1000,pidPlot,'edgecolor','none');
            view(2);
            colormap(s2,cscale_hcr);
            cb=colorbar;
            cb.Ticks=1:9;
            cb.TickLabels=units_str_hcr;
            ylabel('Altitude (km)');
            title(['HCR particle ID']);
            caxis([.5 9.5]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            
            grid on
            box on
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_Flight',num2str(aa),'_pidTemperature_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        startPlot=endPlot;
    end
end