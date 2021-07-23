% Check if temperature correction in the data has been performed

clear all;
close all;

project='spicule';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc1/' project '/checkQC/'];

if strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/moments/100hz/'; %socrates
    highResTempDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
elseif strcmp(project,'cset')
    indir='/scr/rain1/rsfdata/projects/cset/hcr/qc2/cfradial/moments/100hz/'; %cset raw
    highResTempDir='/scr/snow2/rsfdata/projects/cset/hcr/txt/';
elseif strcmp(project,'aristo')
    indir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/'; %aristo raw
    highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
elseif strcmp(project,'otrec')
    indir='/scr/snow1/rsfdata/projects/otrec/hcr/qc1/cfradial/moments/100hz/'; % otrec
    highResTempDir='/scr/snow1/rsfdata/projects/otrec/hcr/txt/';
elseif strcmp(project,'spicule')
    indir='/scr/sleet2/rsfdata/projects/spicule/hcr/cfradial/moments/100hz/'; % spicule field
    highResTempDir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc0/txt/';
else
    disp('Project name not valid.')
    return
end

filedir='~/git/HCR_configuration/projDir/qc/dataProcessing/nsCal/inFiles/';
infile=['cal_' project '.dat'];

inlist=readtable([filedir infile]);
%Infile has start and end dates. In last column it has 0, 1, or 2.
% 0: Use low res temperature data directly from cfradial file, use for
% calculating pod temperature slope
% 1: Use high resolution temperature data when LNA heater is working well,
% i.e. it cycles in regular waves. LNA temperature data will be adjusted for
% time delay between LNA temperatues and power waves. Data used to
% calculate LNA temperature slope and pod temperature slope
% 2: LNA heaters are not working properly. High resolution temperature data
% will be used but LNA temperatures will not be adjusted for time lag. Data
% not used for LNA temperature slope but used for pod temperature
% correction
% 3: Data from lab bench top noise cal. Data will be plotted but not used
% to calculate anything

T0=290;
K=1.38e-23;
ENR_Quinstar = 20.84; %dB

%% Loop through cases

powerDiff=[];
TempAll=[];
dbmvcAll=[];

f1=figure;
set(f1,'Position',[200 100 1200 1000]);
hold on

f2=nan;

if any(inlist{:,13}==1) | any(inlist{:,13}==2)
    legendCell={};
end

colorAll=jet(size(inlist,1));
plotSym='^';

for ii=1:size(inlist,1)
    
    disp(['Case ' num2str(ii) ' from ' num2str(size(inlist,1))]);
    startTime=datetime(inlist{ii,1:6});
    endTime=datetime(inlist{ii,7:12});
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get all relevant data and remove times that are not wanted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dBmVC=[];
    dBZ=[];
    time=[];
    Temp=[];
    timeTemp=[];
    pulseWidth=[];
    
    for jj=1:size(fileList,2)
        disp(['File ',num2str(jj),' from ',num2str(size(fileList,2))]);
        dBmVC=cat(2,dBmVC,ncread(fileList{jj},'DBMVC'));
        dBZ=cat(2,dBZ,ncread(fileList{jj},'DBZ'));
        startTimeIn=ncread(fileList{jj},'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeRead=ncread(fileList{jj},'time')';
        
        timeIn=startTimeFile+seconds(timeRead);
        time=cat(1,time,timeIn');
        if inlist{ii,13}==0 | inlist{ii,13}==3 %If low res temperature is required, read it from the cfradial file
            Temp=[Temp varFromCfRadialString(fileList{jj},'EikTemp')];
            timeTemp=[timeTemp startTimeFile];
        end
        pulseWidthIn=ncread(fileList{jj},'pulse_width');
        pulseWidthChangeIn=find(pulseWidthIn ~= pulseWidthIn(1));
        if ~isempty(pulseWidthChangeIn) % Check if pulse width changes within the file
            disp('Pulse width is not constant!!!');
            return
        else
            pulseWidth1=pulseWidthIn(1);
        end
        pulseWidth=cat(2,pulseWidth,pulseWidth1);
    end
    
    if isempty(dBmVC)
        continue
    end
    % Check if pulse width changes over all files
    pulseWidthChange=find(pulseWidth ~= pulseWidth(1));
    
    if ~isempty(pulseWidthChange)
        disp('Pulse width is not constant!!!');
        return
    end
    
    % Get rid of non noise source cal data
    around90=zeros(size(dBmVC));
    around90(dBmVC>-95 & dBmVC<-85)=1;

    sum90=nansum(around90,1);
    tooSmall=find(sum90<(size(dBmVC,1)*0.9));
    
    dBmVC(:,tooSmall)=[];
    dBZ(:,tooSmall)=[];
    time(tooSmall)=[];
    
    % Get rid of data that is not within the wanted time interval
    outOfBoundsInds=find(time<startTime | time>endTime);
    
    dBmVC(:,outOfBoundsInds)=[];
    dBZ(:,outOfBoundsInds)=[];
    dBmVCmean=nanmean(dBmVC,1);
    dBZmean=nanmean(dBZ,1);
    
    time(outOfBoundsInds)=[];
    
    if isempty(dBmVC)
        continue
    end
    
    %Get high resolution temperature data if wanted
    if inlist{ii,13}==1 | inlist{ii,13}==2
        if strcmp(project,'socrates')
            tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);
            indata=txtTable2matTable(tempFile,',');
        elseif strcmp(project,'cset') | strcmp(project,'aristo') | strcmp(project,'otrec') | strcmp(project,'spicule')
            if ii==1
                if strcmp(project,'cset')
                    tempFile=[highResTempDir,'CSET.temperatures.txt'];
                elseif strcmp(project,'otrec')
                    tempFile=[highResTempDir,'OTREC.temperatures.txt'];
                elseif strcmp(project,'spicule')
                    tempFile=[highResTempDir,'SPICULE.temperatures.txt'];
                elseif strcmp(project,'aristo')
                    tempFile=[highResTempDir,project,'_temps.txt'];
                end
                tempnames={'count','year','month','day','hour','min','sec','unix_time',...
                    'unix_day','XmitterTemp','PloTemp','EikTemp','VLnaTemp','HLnaTemp',...
                    'PolarizationSwitchTemp','RfDetectorTemp','NoiseSourceTemp','Ps28VTemp',...
                    'RdsInDuctTemp','RotationMotorTemp','TiltMotorTemp','CmigitsTemp',...
                    'TailconeTemp','PentekFpgaTemp','PentekBoardTemp'};
                indata=readtable(tempFile);
                indata.Properties.VariableNames=tempnames;
            end
        end
        
        % Remove spaces in variable names if necessary
        varNames=indata.Properties.VariableNames;
        newNames={};
        for aa=1:length(varNames)
            newNames{end+1}=erase(varNames{aa}," ");
        end
        indata.Properties.VariableNames=newNames;
        
        EikTemp=indata.EikTemp;
        PolSwitchTemp=indata.PolarizationSwitchTemp;
        RfDetTemp=indata.RfDetectorTemp;
        NoisSourceTemp=indata.NoiseSourceTemp;
        
        allTemps=cat(2,EikTemp,PolSwitchTemp,RfDetTemp,NoisSourceTemp); % Calculate mean over all relevant temperatures
        Temp=mean(allTemps,2);
        
        %Remove data with the wrong times
        timeTemp=datetime(indata.year,indata.month,indata.day,indata.hour,indata.min,indata.sec);
        outOfTimeInds=find(timeTemp<startTime | timeTemp>endTime);
        Temp(outOfTimeInds)=[];
        timeTemp(outOfTimeInds)=[];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make correlation between mean temperature and powers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dnTime=datenum(time);
    dnTimeTemp=datenum(timeTemp);
    
    % Create time vector where we want to grab the interpolated data
    startWholeMin=dateshift(time(1), 'start', 'minute','nearest');
    endWholeMin=dateshift(time(end), 'start', 'minute','nearest');
    newTimes = startWholeMin:seconds(30):endWholeMin;
    DNnewTimes=datenum(newTimes);
    
    %Interpolate DBMVC data
    nanInd=find(isnan(dnTime));
    dBmVCmean(nanInd)=[];
    dnTime(nanInd)=[];
    time(nanInd)=[];
    
    nanInd=find(isnan(dBmVCmean));
    dBmVCmean(nanInd)=[];
    dnTime(nanInd)=[];
    time(nanInd)=[];
    
    [fitDBMVC,~,mu1]=polyfit(dnTime,dBmVCmean',2);
    newDBMVC=polyval(fitDBMVC,DNnewTimes,[],mu1);
    
    %Interpolate Temperature data
    nanIndT=find(isnan(Temp) | isnan(dnTimeTemp));
    Temp(nanIndT)=[];
    dnTimeTemp(nanIndT)=[];
    
    [fitTemp,~,mu2]=polyfit(dnTimeTemp,Temp,1);
    newTemp=polyval(fitTemp,DNnewTimes,[],mu2);
    
    %Calculate theoretical DBMVC value
    Bn=1/pulseWidth(1)*0.7;
    
    ENR_Corr=10^(ENR_Quinstar/10) + (T0-(newTemp + 273.15))./T0;
    
    Ts_ENR_Corr = ENR_Corr.*T0+T0;
    Pn_ENR_Corr = K.*Ts_ENR_Corr.*Bn;
    
    Pn_ENR_Corr_dB=10.*log10(Pn_ENR_Corr);
    
    powerDiffFile=newDBMVC-(Pn_ENR_Corr_dB+30);
    
    if inlist{ii,13}==0 | inlist{ii,13}==1
        powerDiff=cat(2,powerDiff,powerDiffFile);
        TempAll=cat(2,TempAll,newTemp);
        dbmvcAll=cat(2,dbmvcAll,newDBMVC);
    end
    
    if  ishghandle(f2)
        close(f2);
    end
    
    f2=figure;
    set(f2,'Position',[1500 500 1500 1000]);
    
    subplot(2,1,1)
    hold on;
    
    plot(time,dBmVCmean,'c');
    xlabel('Time [UTC]');
    ylabel('DBMVC [dB]');
    xlim([time(1) time(end)]);
    
    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(time(1),formatOut);
    title([timestring],'interpreter','none');
    
    legend('DBMVC');
    
    subplot(2,1,2)
    hold on;
    
    plot(time,dBZmean,'b');
    plot(time,movmean(dBZmean,600),'g');
    xlabel('Time [UTC]');
    ylabel('DBZ [dBZ]');
    xlim([time(1) time(end)]);
    
    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(time(1),formatOut);
    title([timestring],'interpreter','none');
    
    legend('DBZ','Mean DBZ');
    
    set(f2,'PaperPositionMode','auto')
    print(f2, [figdir 'nsCal_' project '_',timestring],'-dpng','-r0')
    
    if inlist{ii,13}==1 % If LNA heaters are working, adjust time phase shift
        plotSym='o';
    elseif inlist{ii,13}==2
        plotSym='+';
    end
    
    figure(f1);
    scatter(newTemp,powerDiffFile,[],colorAll(ii,:),plotSym);
    drawnow;
    
    legendCell{end+1}=timestring;
end
%% Scatter plot

figure(f1)

xlabel('Temperature [C]');
ylabel('Power difference [dB]');
title([project ' noise source calibration']);

legend(legendCell,'interpreter','none','location','eastoutside');

set(f1,'PaperPositionMode','auto')
print(f1,[figdir 'nsCal_' project],'-dpng','-r0')
