% find minimum reflectivity values
clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
freqData='10hz';

figdir=['/h/eol/romatsch/hcrCalib/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freqData);

altAll=[];
hAll=[];
vAll=[];
%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6))+minutes(30);
    endTimeAll=datetime(caseList(ii,7:12))-minutes(30);
    
    altFlight=[];
    timeFlight=[];
    elevFlight=[];
    hFlight=[];
    vFlight=[];

    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
                
        % Desired variables. The variable name comies after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        data.DBMVC=[];
        data.DBMHX=[];
        data.DBZ=[];
        
        dataVars=fieldnames(data);
              
        %% Load data
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            startTime=endTime;
            continue
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        if isempty(data.time)
            disp('No data found.');
            startTime=endTime;
            continue
        end
        
        % Check if all variables were found
        for kk=1:length(dataVars)
            if ~isfield(data,dataVars{kk})
                dataVars{kk}=[];
            end
        end
                                 
        % Remove up pointing
        data.DBZ(:,data.elevation>-89.98)=nan;
        data.elevation(:,data.elevation>-89.98)=nan;
        
        % Surface ind
        [linInd rowInd rangeToSurf] = hcrSurfInds(data);
        hsurf=data.DBMHX(linInd);
        vsurf=data.DBMVC(linInd);
        
        hsurf(hsurf<-90)=nan;
        vsurf(vsurf<-60)=nan;
                
        altFlight=cat(2,altFlight,data.altitude);
        timeFlight=cat(2,timeFlight,data.time);
        elevFlight=cat(2,elevFlight,data.elevation);
        
        hFlight=cat(2,hFlight,hsurf);
        vFlight=cat(2,vFlight,vsurf);
        
        startTime=endTime;
        
    end
    
    altAll=cat(2,altAll,altFlight);
    hAll=cat(2,hAll,hFlight);
    vAll=cat(2,vAll,vFlight);
    
    %% Noise figures
    
    close all
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 1300 1200]);
    
    subplot(3,1,1)
    hold on
    plot(timeFlight,vFlight,'linewidth',1.5);
    plot(timeFlight,hFlight,'linewidth',1.5);
    
    ylabel('Power (dB)');
    xlim([timeFlight(1) timeFlight(end)]);
    title([project,' flight ',num2str(ii)]);
    
    legend('Surf DBMVC','Surf DBMHX');
    grid on
    
    subplot(3,1,2)
    hold on
    plot(timeFlight,vFlight-hFlight,'linewidth',1.5);
    
    ylabel('Power (dB)');
    xlim([timeFlight(1) timeFlight(end)]);
    grid on
    
    ylim([15 35]);
    legend('Surf DBMVC - Surf DBMHX');
    
    subplot(3,1,3)
    hold on
    plot(timeFlight,altFlight./1000,'linewidth',1.5);
    
    ylabel('Altitude (km)');
    xlim([timeFlight(1) timeFlight(end)]);
    grid on
    
    legend('Altitude');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'extremePower/',project,'/','timeSeriesVH_',project,'_RF',num2str(ii)],'-dpng','-r0');  
    
    % Scatter
    f3=figure('DefaultAxesFontSize',12);
    set(f3,'Position',[200 500 800 800]);
          
    scatter(vFlight-hFlight,altFlight./1000);
    grid on
    
    title([project,' flight ',num2str(ii)]);
    xlabel('Surf DBMVC - Surf DBMHX');
    ylabel('Altitude (km)');
    xlim([0 40]);
    ylim([0 15]);
                        
     set(gcf,'PaperPositionMode','auto')
    print(f3,[figdir,'extremePower/',project,'/','scatterVH_',project,'_RF',num2str(ii)],'-dpng','-r0');  
    
    
end

close all

f4=figure('DefaultAxesFontSize',12);
set(f4,'Position',[200 500 800 800]);

scatter(vAll-hAll,altAll./1000);
grid on

title([project]);
xlabel('Surf DBMVC - Surf DBMHX');
ylabel('Altitude (km)');
xlim([0 40]);
ylim([0 15]);

set(gcf,'PaperPositionMode','auto')
print(f4,[figdir,'extremePower/',project,'/','scatterVH_',project],'-dpng','-r0');