% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc0';
freqData='10hz_combined';
qcVersion='';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'sensitivity/'];

edges=-100:0.5:100;
edgesS=-30:0.5:100;

dbzAllS=zeros(1,length(edges)-1);
dbz1kmS=zeros(1,length(edges)-1);
dbz5kmS=zeros(1,length(edges)-1);

snrAllS=zeros(1,length(edgesS)-1);
snr1kmS=zeros(1,length(edgesS)-1);
snr5kmS=zeros(1,length(edgesS)-1);

dbzAllL=zeros(1,length(edges)-1);
dbz1kmL=zeros(1,length(edges)-1);
dbz5kmL=zeros(1,length(edges)-1);

snrAllL=zeros(1,length(edgesS)-1);
snr1kmL=zeros(1,length(edgesS)-1);
snr5kmL=zeros(1,length(edgesS)-1);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
            
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+hours(1);

        data=[];
                   
        data.DBZ_short=[];
        data.DBZ_long=[];
        data.SNRVC_short=[];
        data.SNRVC_long=[];
                              
        %% Load data
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
                
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
                
        asl=HCRrange2asl(data.range,data.elevation,data.altitude);
        
        dbzStemp=data.DBZ_short;
        dbzLtemp=data.DBZ_long;
        snrStemp=data.SNRVC_short;
        snrLtemp=data.SNRVC_long;
                
        % Remove bang
        dbzStemp(1:14,:)=nan;
        snrStemp(1:14,:)=nan;
        dbzLtemp(1:14,:)=nan;
        snrLtemp(1:14,:)=nan;
        % Remove data that is 300 m or below from the surface
        dbzStemp(asl<300)=nan;
        snrStemp(asl<300)=nan;
        dbzLtemp(asl<300)=nan;
        snrLtemp(asl<300)=nan;
                
        histDBZs=histcounts(dbzStemp,edges);
        dbzAllS=dbzAllS+histDBZs;
        histDBZl=histcounts(dbzLtemp,edges);
        dbzAllL=dbzAllL+histDBZl;
        
        histSNRs=histcounts(snrStemp,edgesS);
        snrAllS=snrAllS+histSNRs;
        histSNRl=histcounts(snrLtemp,edgesS);
        snrAllL=snrAllL+histSNRl;
        
        range1km=find(data.range>=990 & data.range <=1010);
        range5km=find(data.range>=4990 & data.range <=5010);
        
        histDBZ1kmS=histcounts(dbzStemp(range1km),edges);
        histDBZ5kmS=histcounts(dbzStemp(range5km),edges);
        histDBZ1kmL=histcounts(dbzLtemp(range1km),edges);
        histDBZ5kmL=histcounts(dbzLtemp(range5km),edges);
        
        histSNR1kmS=histcounts(snrStemp(range1km),edgesS);
        histSNR5kmS=histcounts(snrStemp(range5km),edgesS);
        histSNR1kmL=histcounts(snrLtemp(range1km),edgesS);
        histSNR5kmL=histcounts(snrLtemp(range5km),edgesS);
        
        dbz1kmS=dbz1kmS+histDBZ1kmS;
        dbz5kmS=dbz5kmS+histDBZ5kmS;
        dbz1kmL=dbz1kmL+histDBZ1kmL;
        dbz5kmL=dbz5kmL+histDBZ5kmL;
        
        snr1kmS=snr1kmS+histSNR1kmS;
        snr5kmS=snr5kmS+histSNR5kmS;
        snr1kmL=snr1kmL+histSNR1kmL;
        snr5kmL=snr5kmL+histSNR5kmL;
        
        startTime=endTime;
    end
end
    %% All
close all
    f1 = figure('Position',[200 500 1000 1250],'DefaultAxesFontSize',12);

t = tiledlayout(4,2,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);
    
    bar(edges(1:end-1)+0.25,dbzAllS,1);
    xlabel('dBZ');
    title([project,' DBZ short']);
    xlim([-65 20]);
    grid on
    box on
        
    s2=nexttile(2);
    bar(edges(1:end-1)+0.25,dbzAllS,1);
    xlabel('dBZ');
    title(['DBZ short lower end']);
    xlim([-65 -40]);
    grid on
    box on

    s3=nexttile(3);
    
    bar(edges(1:end-1)+0.25,dbzAllL,1);
    xlabel('dBZ');
    title([project,' DBZ long']);
    xlim([-65 20]);
    grid on
    box on
        
    s4=nexttile(4);
    bar(edges(1:end-1)+0.25,dbzAllL,1);
    xlabel('dBZ');
    title(['DBZ long lower end']);
    xlim([-65 -40]);
    grid on
    box on
    
    s5=nexttile(5);
    bar(edgesS(1:end-1)+0.25,snrAllS,1);
    xlabel('dB');
    title(['SNRVC short']);
    xlim([-25 70]);
    grid on
    box on
    
    s6=nexttile(6);
    bar(edgesS(1:end-1)+0.25,snrAllS,1);
    xlabel('dB');
    title(['SNRVC short lower end']);
    xlim([-25 -5]);
    grid on
    box on

    s7=nexttile(7);
    bar(edgesS(1:end-1)+0.25,snrAllL,1);
    xlabel('dB');
    title(['SNRVC long']);
    xlim([-25 70]);
    grid on
    box on
    
    s8=nexttile(8);
    bar(edgesS(1:end-1)+0.25,snrAllL,1);
    xlabel('dB');
    title(['SNRVC long lower end']);
    xlim([-25 -5]);
    grid on
    box on
            
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_minDBZ'],'-dpng','-r0')
    
    %% 1km
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 1300 800]);
       
    subplot(2,2,1)
    
    bar(edges(1:end-1)+0.25,dbz1kmS,1);
    xlabel('dBZ');
    title([project,' DBZ at 1 km range']);
    xlim([-50 30]);
    
    subplot(2,2,2)
    
    bar(edges(1:end-1)+0.25,dbz1kmS,1);
    xlabel('dBZ');
    title(['DBZ lower end']);
    xlim([-50 -20]);
    
    subplot(2,2,3)
    
    bar(edgesS(1:end-1)+0.25,snr1kmS,1);
    xlabel('dB');
    title(['SNRVC at 1 km range']);
    xlim([-25 50]);
    
    subplot(2,2,4)
    
    bar(edgesS(1:end-1)+0.25,snr1kmS,1);
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
    
    bar(edges(1:end-1)+0.25,dbz5kmS,1);
    xlabel('dBZ');
    title([project,' DBZ at 5 km range']);
    xlim([-40 40]);
    
    subplot(2,2,2)
    
    bar(edges(1:end-1)+0.25,dbz5kmS,1);
    xlabel('dBZ');
    title(['DBZ lower end']);
    xlim([-40 -10]);
    
    subplot(2,2,3)
    
    bar(edgesS(1:end-1)+0.25,snr5kmS,1);
    xlabel('dB');
    title(['SNRVC at 5 km range']);
    xlim([-25 50]);
    
    subplot(2,2,4)
    
    bar(edgesS(1:end-1)+0.25,snr5kmS,1);
    xlabel('dB');
    title(['SNRVC lower end']);
    xlim([-25 -5]);
            
    set(gcf,'PaperPositionMode','auto')
    print(f3,[figdir,project,'_minDBZ_5km'],'-dpng','-r0')
    