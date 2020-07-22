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

edges=0:0.1:15;

histAll=[];
%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
    
    altFlight=zeros(1,length(edges)-1);

    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
                
        % Desired variables. The variable name comies after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        data.dummy=[];
                
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
        data.altitude(data.elevation>0)=nan;
                 
        altFlight=altFlight+histcounts(data.altitude./1000,edges);
        
        startTime=endTime;
        
    end
    histAll=cat(1,histAll,altFlight);    
end

%% Plot flights
for ii=1:size(histAll,1)
    
    close all
    
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 800 800]);
    
    barh(edges(1:end-1),histAll(ii,:)./10,1);
    grid on
    set(gca, 'XScale', 'log')
        
    title([project,' flight ',num2str(ii)]);
    xlabel('Seconds');
    ylabel('Altitude (km)');
    xlim([1 10e+4]);
    ylim([0 14]);
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'extremePower/',project,'/','timeVSalt_',project,'_flight',num2str(ii)],'-dpng','-r0');
end

f1=figure('DefaultAxesFontSize',12);
set(f1,'Position',[200 500 800 800]);

barh(edges(1:end-1),sum(histAll,1)./10,1);
grid on
set(gca, 'XScale', 'log')

title([project]);
xlabel('Seconds');
ylabel('Altitude (km)');
xlim([1 10e+5]);
ylim([0 14]);

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,'extremePower/',project,'/','timeVSalt_',project],'-dpng','-r0');