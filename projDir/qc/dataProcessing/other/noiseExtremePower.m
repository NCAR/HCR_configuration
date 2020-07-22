% find minimum reflectivity values
clear all;
close all;

addpath(genpath('/h/eol/romatsch/gitPriv/utils/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
freqData='100hz';

figdir=['/h/eol/romatsch/hcrCalib/'];
formatOut = 'yyyymmdd_HHMM';

infile=['/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freqData);

edgesNoise=-110:0.1:-90;
edgesMax=-120:1:120;

vAll=zeros(1,length(edgesNoise)-1);
hAll=zeros(1,length(edgesNoise)-1);

noiseMedV=[];
noiseMedH=[];
noiseSumV=0;
noiseSumH=0;
noiseSumV2=0;
noiseSumH2=0;
time10minV=[];
time10minH=[];

flightMeanVall=[];
flightMeanHall=[];
flightStdVall=[];
flightStdHall=[];

NV=0;
NH=0;

maxV=zeros(1,length(edgesMax)-1);
maxH=zeros(1,length(edgesMax)-1);
maxHatV=zeros(1,length(edgesMax)-1);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    startTime=datetime(caseList(ii,1:6));
    endTimeAll=datetime(caseList(ii,7:12));
    
    vFlight=zeros(1,length(edgesNoise)-1);
    hFlight=zeros(1,length(edgesNoise)-1);
    
    noiseMedVF=[];
    noiseMedHF=[];
    noiseSumVF=0;
    noiseSumHF=0;
    noiseSumV2F=0;
    noiseSumH2F=0;
    
    NVF=0;
    NHF=0;
    
    maxVflight=[];
    maxHflight=[];
    maxHatVflight=[];
    
    altFlight=[];
    timeFlight=[];
    
    vc1dF=[];
    hx1dF=[];
    
    % Go through hours
    while startTime<endTimeAll
        disp(startTime);
        
        endTime=startTime+minutes(10);
                
        % Desired variables. The variable name comies after the . and must be spelled exactly
        % like in the CfRadial file
        if exist('data')
            clear data
        end
        
        data.DBMVC=[];
        data.DBMHX=[];
               
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
                    
        altFlight=cat(2,altFlight,data.altitude);
        timeFlight=cat(2,timeFlight,data.time);
        
        %% Calculate noise
        
        vctemp=data.DBMVC;
        hxtemp=data.DBMHX;
        
        % Remove down pointing
        vctemp(:,data.elevation<0)=nan;
        hxtemp(:,data.elevation<0)=nan;
        
        % Find clouds
        cloud=data.DBMVC;
        cloud(1:12,:)=nan;
        cloud(cloud<-100)=nan;
        
        cloudMask=zeros(size(data.DBMVC));
        cloudMask(~isnan(cloud))=1;
        
        cloudSum=sum(cloudMask,1);
        meanCloud=movmean(cloudSum,100);
        
        vctemp(:,meanCloud>4)=nan;
        hxtemp(:,meanCloud>4)=nan;
        vctemp(:,cloudSum>6)=nan;
        hxtemp(:,cloudSum>6)=nan;
        
        vc100=vctemp(671:770,:);
        hx100=hxtemp(671:770,:);
        
        vc100lin=10.^(vc100./10);
        hx100lin=10.^(hx100./10);
        
        meanVC=nanmean(vc100lin,1);
        meanHX=nanmean(hx100lin,1);
        
        vc1d=log10(meanVC).*10;
        hx1d=log10(meanHX).*10;
        
        vc1dF=cat(2,vc1dF,vc1d);
        hx1dF=cat(2,hx1dF,hx1d);
        
        histV=histcounts(vc100,edgesNoise);
        vFlight=vFlight+histV;
        vAll=vAll+histV;
        histH=histcounts(hx100,edgesNoise);
        hFlight=hFlight+histH;
        hAll=hAll+histH;
        
        vc100vec=reshape(vc100,1,[]);
        hx100vec=reshape(hx100,1,[]);
        
        vc100vec(isnan(vc100vec))=[];
        hx100vec(isnan(hx100vec))=[];
        
        if ~isempty(vc100vec)
            noiseMedVF=[noiseMedVF,mean(vc100vec)];
            noiseSumV2F=noiseSumV2F+sum(vc100vec.^2);
            noiseSumVF=noiseSumVF+sum(vc100vec);
            NVF=NVF+length(vc100vec);
            time10minV=[time10minV startTime+minutes(5)];
        end
        if ~isempty(hx100vec)
            noiseMedHF=[noiseMedHF,mean(hx100vec)];
            noiseSumH2F=noiseSumH2F+sum(hx100vec.^2);
            noiseSumHF=noiseSumHF+sum(hx100vec);
            NHF=NHF+length(hx100vec);
            time10minH=[time10minH startTime+minutes(5)];
        end
        
        %% Extreme power
        
        vctemp=data.DBMVC;
        hxtemp=data.DBMHX;
        
        % Remove bang
        vctemp(1:15,:)=nan;
        hxtemp(1:15,:)=nan;
        
        % Remove up pointing
        vctemp(:,data.elevation>0)=nan;
        hxtemp(:,data.elevation>0)=nan;
        
        [maxVf maxVI]=nanmax(vctemp,[],1);
        maxVflight=cat(2,maxVflight,maxVf);
        
        maxHf=nanmax(hxtemp,[],1);
        maxHflight=cat(2,maxHflight,maxHf);
        
        linV=sub2ind(size(hxtemp),maxVI,(1:1:size(hxtemp,2)));
        maxHatVflight=cat(2,maxHatVflight,hxtemp(linV));
        
        startTime=endTime;
                
    end
    
    noiseMedV=cat(2,noiseMedV,noiseMedVF);
    noiseMedH=cat(2,noiseMedH,noiseMedHF);
    noiseSumV=noiseSumV+noiseSumVF;
    noiseSumH=noiseSumH+noiseSumHF;
    noiseSumV2=noiseSumV2+noiseSumV2F;
    noiseSumH2=noiseSumH2+noiseSumH2F;
    
    NV=NV+NVF;
    NH=NH+NHF;
    
    flightMeanV=mean(noiseMedVF);
    flightMeanH=mean(noiseMedHF);
    
    flightStdV=sqrt(noiseSumV2F/NVF-(noiseSumVF/NVF)^2);
    flightStdH=sqrt(noiseSumH2F/NHF-(noiseSumHF/NHF)^2);
    
    flightMeanVall=[flightMeanVall,flightMeanV];
    flightMeanHall=[flightMeanHall,flightMeanH];
    flightStdVall=[flightStdVall,flightStdV];
    flightStdHall=[flightStdHall,flightStdH];
    
    %% Noise figures
    
    close all
    f1=figure('DefaultAxesFontSize',12);
    set(f1,'Position',[200 500 1300 800]);
       
    subplot(2,1,1)
    
    bar(edgesNoise(1:end-1),vFlight,1);
    xlabel('DBMVC noise');
    xlim([-106 -98]);
    title([project,' flight ',num2str(ii)]);
    grid on
    
    text(-105.9,max(vFlight)/10,['Mean: ',num2str(flightMeanV),', Std: ',num2str(flightStdV)],'fontsize',12);
    
    subplot(2,1,2)
    
    bar(edgesNoise(1:end-1),hFlight,1);
    xlabel('DBMHX noise');
    xlim([-106 -98]);
    grid on
    
    text(-105.9,max(hFlight)/10,['Mean: ',num2str(flightMeanH),', Std: ',num2str(flightStdH)],'fontsize',12);
                    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,'noise/',project,'/hist_',project,'_RF',num2str(ii)],'-dpng','-r0');  
    
    % Time series
    f3=figure('DefaultAxesFontSize',12);
    set(f3,'Position',[200 500 1800 600]);
          
    hold on
    plot(timeFlight,vc1dF,'linewidth',2);
    plot(timeFlight,hx1dF,'linewidth',2);
    ylabel('Noise (dB)');
    ylim([-106 -98.5])
    yticks(-106:-99);
    
    yyaxis right
    plot(timeFlight,altFlight./1000,'linewidth',2);
    ylabel('Altitune (km)');
    xlim([timeFlight(1) timeFlight(end)]);
    title([project,' flight ',num2str(ii)]);
    ylim([0 15])
    yticks(0:15)
    grid on
                        
    legend('DBMVC noise','DBMHX noise','Altitude')
    set(gcf,'PaperPositionMode','auto')
    print(f3,[figdir,'noise/',project,'/timeSeries_',project,'_RF',num2str(ii)],'-dpng','-r0');   
    
    %% Extreme power
    
    linMaxV=10.^(maxVflight./10);
    meanMaxV=movmean(linMaxV,100);
    meanVdb=log10(meanMaxV).*10;
    
    histVmax=histcounts(meanVdb,edgesMax);
    maxV=maxV+histVmax;
    
    linMaxH=10.^(maxHflight./10);
    meanMaxH=movmean(linMaxH,100);
    meanHdb=log10(meanMaxH).*10;
    
    histHmax=histcounts(meanHdb,edgesMax);
    maxH=maxH+histHmax;
    
    linMaxHatV=10.^(maxHatVflight./10);
    meanMaxHatV=movmean(linMaxHatV,100);
    meanHatVdb=log10(meanMaxHatV).*10;
    
    histHatVmax=histcounts(meanHatVdb,edgesMax);
    maxHatV=maxHatV+histHatVmax;    
    
    %% Hist extreme power
    
    f2=figure('DefaultAxesFontSize',12);
    set(f2,'Position',[200 500 1300 1000]);
       
    subplot(2,2,1)
    
    bar(edgesMax(1:end-1),histVmax,1);
    xlabel('Max DBMVC 1s run mean');
    xlim([-110 0]);
    ylim([0 9e+5]);
    title([project,' flight ',num2str(ii)]);
    
    grid on
    
    subplot(2,2,2)
    
    bar(edgesMax(1:end-1),histVmax,1);
    xlabel('Max DBMVC 1s run mean high end');
    xlim([-50 0]);
    ylim([0 9e+5]);
    set(gca,'YScale','log')
        
    grid on
 
    subplot(2,2,3)
    
    bar(edgesMax(1:end-1)+27,histHatVmax,1);
    xlabel('Max DBMHX at max DBMVC gate 1s run mean + 27');
    xlim([-110 0]);
    ylim([0 9e+5]);
    
    grid on
    
    subplot(2,2,4)
    
    bar(edgesMax(1:end-1)+27,histHatVmax,1);
    xlabel('Max DBMHX at max DBMVC gate 1s run mean + 27');
    xlim([-50 0]);
    ylim([0 9e+5]);
    set(gca,'YScale','log')
    
    grid on
                    
    set(gcf,'PaperPositionMode','auto')
    print(f2,[figdir,'extremePower/',project,'/hist_',project,'_RF',num2str(ii)],'-dpng','-r0');
end
%% 

MeanV=mean(noiseMedV);
MeanH=mean(noiseMedH);
    
StdV=sqrt(noiseSumV2/NV-(noiseSumV/NV)^2);
StdH=sqrt(noiseSumH2/NH-(noiseSumH/NH)^2);
    
%% All

close all
f4=figure('DefaultAxesFontSize',12);
set(f4,'Position',[200 500 1300 800]);

subplot(2,1,1)

bar(edgesNoise(1:end-1),vAll,1);
xlabel('DBMVC noise');
xlim([-106 -98]);
title([project]);
grid on

text(-105.9,max(vAll)/10,['Mean: ',num2str(MeanV),', Std: ',num2str(StdV)],'fontsize',12);

subplot(2,1,2)

bar(edgesNoise(1:end-1),hAll,1);
xlabel('DBMHX noise');
xlim([-106 -98]);
grid on

text(-105.9,max(hAll)/10,['Mean: ',num2str(MeanH),', Std: ',num2str(StdH)],'fontsize',12);

set(gcf,'PaperPositionMode','auto')
print(f4,[figdir,'noise/',project,'/hist_',project],'-dpng','-r0');

f6=figure('DefaultAxesFontSize',12);
set(f6,'Position',[200 500 1300 800]);

subplot(2,1,1)
hold on
plot(time10minV,noiseMedV,'.')
plot(time10minH,noiseMedH,'.')
ylabel('DBM noise');
ylim([-106 -100]);
xlim([time10minV(1),time10minV(end)]);
title([project]);

legend('DBMVC noise','DBMHX noise','location','southeast');

subplot(2,1,2)
hold on
errorbar((1:1:ii),flightMeanVall,flightStdVall,'linewidth',1.5);
errorbar((1:1:ii),flightMeanHall,flightStdHall,'linewidth',1.5);
xlabel('Flight');
ylabel('DBM noise');
ylim([-106 -100]);
xticks(1:1:ii);
xlim([1,ii]);

legend('DBMVC noise','DBMHX noise','location','southeast');

set(gcf,'PaperPositionMode','auto')
print(f6,[figdir,'noise/',project,'/noiseTime_',project],'-dpng','-r0');

%% Hist extreme power

f5=figure('DefaultAxesFontSize',12);
set(f5,'Position',[200 500 1300 1000]);

subplot(2,2,1)

bar(edgesMax(1:end-1),maxV,1);
xlabel('Max DBMVC 1s run mean');
xlim([-110 0]);
ylim([0 6e+6]);
title([project]);
grid on

subplot(2,2,2)

bar(edgesMax(1:end-1),maxV,1);
xlabel('Max DBMVC 1s run mean high end');
xlim([-50 0]);
ylim([0 6e+6]);
set(gca,'YScale','log')
grid on

subplot(2,2,3)

bar(edgesMax(1:end-1)+27,maxHatV,1);
xlabel('Max DBMHX at max DBMVC gate 1s run mean + 27');
xlim([-110 0]);
ylim([0 6e+6]);
grid on

subplot(2,2,4)

bar(edgesMax(1:end-1)+27,maxHatV,1);
xlabel('Max DBMHX at max DBMVC gate 1s run mean + 27');
xlim([-50 0]);
ylim([0 6e+6]);
set(gca,'YScale','log')
grid on

set(gcf,'PaperPositionMode','auto')
print(f5,[figdir,'extremePower/',project,'/hist_',project],'-dpng','-r0');
%% 

edgesMax27=edgesMax(1:end-1)+27;

ind20=find(edgesMax27>-20);
sum20=nansum(maxHatV(ind20))/100;
min20=sum20/60;

ind30=find(edgesMax27>-30);
sum30=nansum(maxHatV(ind30))/100;
min30=sum30/60;

disp(['Over -20 dB: ',num2str(sum20),' seconds, ',num2str(min20),' minutes.']);
disp(['Over -30 dB: ',num2str(sum30),' seconds, ',num2str(min30),' minutes.']);