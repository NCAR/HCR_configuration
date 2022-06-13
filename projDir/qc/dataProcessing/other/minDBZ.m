% find minimum reflectivity values
clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
freqData='10hz';
qcVersion='v1.1';

figdir=['/h/eol/romatsch/hcrCalib/sensitivity/'];
formatOut = 'yyyymmdd_HHMM';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

edges=-100:0.5:100;
edgesS=-30:0.5:100;

dbzAll=zeros(1,length(edges)-1);
dbz1km=zeros(1,length(edges)-1);
dbz5km=zeros(1,length(edges)-1);

snrAll=zeros(1,length(edgesS)-1);
snr1km=zeros(1,length(edgesS)-1);
snr5km=zeros(1,length(edgesS)-1);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
            
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);
                
        % Desired variables. The variable name comies after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        data.DBZ=[];
        if strcmp(quality,'qc2')
            data.SNR=[];
        else
            data.SNRVC=[];
        end
        
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
        
        asl=HCRrange2asl(data.range,data.elevation,data.altitude);
        
        dbztemp=data.DBZ;
        if isfield(data,'SNR')
            snrtemp=data.SNR;
        else
            snrtemp=data.SNRVC;
        end
        
        % Remove bang
        dbztemp(1:14,:)=nan;
        snrtemp(1:14,:)=nan;
        % Remove data that is 300 m or below from the surface
        dbztemp(asl<300)=nan;
        snrtemp(asl<300)=nan;
        
        histDBZ=histcounts(dbztemp,edges);
        dbzAll=dbzAll+histDBZ;
        
        histSNR=histcounts(snrtemp,edgesS);
        snrAll=snrAll+histSNR;
        
        range1km=find(data.range>=990 & data.range <=1010);
        range5km=find(data.range>=4990 & data.range <=5010);
        
        histDBZ1km=histcounts(dbztemp(range1km),edges);
        histDBZ5km=histcounts(dbztemp(range5km),edges);
        
        histSNR1km=histcounts(snrtemp(range1km),edgesS);
        histSNR5km=histcounts(snrtemp(range5km),edgesS);
        
        dbz1km=dbz1km+histDBZ1km;
        dbz5km=dbz5km+histDBZ5km;
        
        snr1km=snr1km+histSNR1km;
        snr5km=snr5km+histSNR5km;
        
        startTime=endTime;
    end
end
    %% All
    close all
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 1300 800]);
       
    subplot(2,2,1)
    
    bar(edges(1:end-1)+0.25,dbzAll,1);
    xlabel('dBZ');
    title([project,' DBZ']);
    xlim([-70 40]);
    
    subplot(2,2,2)
    
    bar(edges(1:end-1)+0.25,dbzAll,1);
    xlabel('dBZ');
    title(['DBZ lower end']);
    xlim([-70 -20]);
    
    subplot(2,2,3)
    
    bar(edgesS(1:end-1)+0.25,snrAll,1);
    xlabel('dB');
    title(['SNRVC']);
    xlim([-25 70]);
    
    subplot(2,2,4)
    
    bar(edgesS(1:end-1)+0.25,snrAll,1);
    xlabel('dB');
    title(['SNRVC lower end']);
    xlim([-25 -5]);
            
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_minDBZ'],'-dpng','-r0')
    
    %% 1km
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 1300 800]);
       
    subplot(2,2,1)
    
    bar(edges(1:end-1)+0.25,dbz1km,1);
    xlabel('dBZ');
    title([project,' DBZ at 1 km range']);
    xlim([-50 30]);
    
    subplot(2,2,2)
    
    bar(edges(1:end-1)+0.25,dbz1km,1);
    xlabel('dBZ');
    title(['DBZ lower end']);
    xlim([-50 -20]);
    
    subplot(2,2,3)
    
    bar(edgesS(1:end-1)+0.25,snr1km,1);
    xlabel('dB');
    title(['SNRVC at 1 km range']);
    xlim([-25 50]);
    
    subplot(2,2,4)
    
    bar(edgesS(1:end-1)+0.25,snr1km,1);
    xlabel('dB');
    title(['SNRVC lower end']);
    xlim([-25 -5]);
            
    set(gcf,'PaperPositionMode','auto')
    print(f2,[figdir,project,'_minDBZ_1km'],'-dpng','-r0')
    savefig([figdir,project,'_minDBZ_1km.fig'])
    
    %% 5km
    
    f3=figure('DefaultAxesFontSize',12);
    set(f3,'Position',[200 500 1300 800]);
       
    subplot(2,2,1)
    
    bar(edges(1:end-1)+0.25,dbz5km,1);
    xlabel('dBZ');
    title([project,' DBZ at 5 km range']);
    xlim([-40 40]);
    
    subplot(2,2,2)
    
    bar(edges(1:end-1)+0.25,dbz5km,1);
    xlabel('dBZ');
    title(['DBZ lower end']);
    xlim([-40 -10]);
    
    subplot(2,2,3)
    
    bar(edgesS(1:end-1)+0.25,snr5km,1);
    xlabel('dB');
    title(['SNRVC at 5 km range']);
    xlim([-25 50]);
    
    subplot(2,2,4)
    
    bar(edgesS(1:end-1)+0.25,snr5km,1);
    xlabel('dB');
    title(['SNRVC lower end']);
    xlim([-25 -5]);
            
    set(gcf,'PaperPositionMode','auto')
    print(f3,[figdir,project,'_minDBZ_5km'],'-dpng','-r0')
    