% Perform noise cal

clear all;
close all;

dispFig='off';

project='spicule';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc1/' project '/lnaTemps/'];

if strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; %socrates
    highResTempDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
elseif strcmp(project,'cset')
    indir='/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/'; %cset raw
    highResTempDir='/scr/snow2/rsfdata/projects/cset/hcr/txt/';
elseif strcmp(project,'aristo')
    indir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/'; %aristo raw
    highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
elseif strcmp(project,'otrec')
    indir='/scr/snow1/rsfdata/projects/otrec/hcr/qc0/cfradial/moments/100hz/'; % otrec field
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

f1=figure('visible',dispFig);
set(f1,'Position',[200 100 1200 1000]);
hold on

f2=nan;
f3=nan;
f4=nan;

if any(inlist{:,13}==1) | any(inlist{:,13}==2)
    f5=figure('visible',dispFig);
    set(f5,'Position',[200 100 1200 1000]);
    hold on
    legendCell={};
    
    timeLagAll=[];
    powTempSlope=[];
    lnaTempsAll=[];
end

f6=figure('visible',dispFig);
set(f6,'Position',[200 100 1200 1000]);
hold on
legendCell={};

colorAll=jet(size(inlist,1));
plotSym='^';

for ii=1:size(inlist,1)
    
    %% Load HCR data
    
    disp(['Case ' num2str(ii) ' from ' num2str(size(inlist,1))]);
    startTime=datetime(inlist{ii,1:6});
    endTime=datetime(inlist{ii,7:12});
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get all relevant data and remove times that are not wanted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dBmVC=[];
    dBZ=[];
    range=[];
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
        
    %% Load High res temperature data
    
    disp('Loading temperature data ...')
    %Get high resolution temperature data if wanted
    if inlist{ii,13}==1 | inlist{ii,13}==2
       if strcmp(project,'socrates')
           tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);
           indata=txtTable2matTable(tempFile,',');
       elseif strcmp(project,'cset') | strcmp(project,'aristo') | strcmp(project,'otrec') | strcmp(project,'spicule')
           if ~exist('indata')
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
       VLnaTemp=indata.VLnaTemp;
       HLnaTemp=indata.HLnaTemp;
       
       allTemps=cat(2,EikTemp,PolSwitchTemp,RfDetTemp,NoisSourceTemp); % Calculate mean over all relevant temperatures
       Temp=mean(allTemps,2);
       
       %Remove data with the wrong times
       timeTemp=datetime(indata.year,indata.month,indata.day,indata.hour,indata.min,indata.sec);
       outOfTimeInds=find(timeTemp<startTime | timeTemp>endTime);
       Temp(outOfTimeInds)=[];
       VLnaTemp(outOfTimeInds)=[];
       HLnaTemp(outOfTimeInds)=[];
       timeTemp(outOfTimeInds)=[];
    end
    
    %% Make correlation between mean temperature and powers
            
    disp('Finding temperature-power correlation ...');
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
    timeTemp(nanIndT)=[];
    VLnaTemp(nanIndT)=[];
    HLnaTemp(nanIndT)=[];
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
    if  ishghandle(f3)
        close(f3);
    end
    if  ishghandle(f4)
        close(f4);
    end
    
    % Debug figure
%     figure
%     plot(time,dBZmean);
%     xlim([time(1) time(end)]);
%     hold on
%     plot(time,movmean(dBZmean,600));
    % %%%%%%%%%%%
    
    f2=figure('visible',dispFig);
    set(f2,'Position',[1500 500 1500 500]);
    
    hold on;
    
    plot(time,dBmVCmean,'c');
    plot(newTimes,newDBMVC,'-bo','MarkerFaceColor','b');
    plot(newTimes,Pn_ENR_Corr_dB+30,'-ro','MarkerFaceColor','r');
    xlabel('Time [UTC]');
    ylabel('DBMVC [dB]');
    xlim([time(1) time(end)]);
    
    yyaxis right
    plot(timeTemp,Temp,'-go','MarkerFaceColor','g')
    plot(newTimes,newTemp,'-ko','MarkerFaceColor','k')
    ylabel('Temperature [C]');
    
    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(time(1),formatOut);
    title([timestring],'interpreter','none');
    
    legend('DBMVC raw','DBMVC interp','Noise Source','Temperature raw','Temperature interp');
    
    
    set(f2,'PaperPositionMode','auto')
    print(f2, [figdir 'nsCal_' project '_',timestring],'-dpng','-r0')
    
    %% Plot LNA temperature vs powers for high res temp data
    
    if  inlist{ii,13}==1 | inlist{ii,13}==2     
       
        % Smooth LNA temperature data
        meanVLNAt=movmean(VLnaTemp,20);
        powers=timetable(time,dBmVCmean');
        vlnatemps=timetable(timeTemp,meanVLNAt);

        synchData=synchronize(powers,vlnatemps,'secondly','mean');
        
        synchData=rmmissing(synchData);
                
        if inlist{ii,13}==1 % If LNA heaters are working, adjust time phase shift
            %caculate phase shift
            pr=max(synchData.Var1)-min(synchData.Var1);
            tr=max(synchData.meanVLNAt)-min(synchData.meanVLNAt);
            
            pProm=max([pr*0.2 0.05]);
            tProm=max([tr*0.2 0.3]);
            
            [peaksP locsP]=findpeaks(synchData.Var1,'MinPeakDistance',30,'MinPeakProminence',pProm);
            [peaksT locsT]=findpeaks(synchData.meanVLNAt,'MinPeakDistance',30,'MinPeakProminence',tProm);
            
            [valleysP locsP2]=findpeaks(-synchData.Var1,'MinPeakDistance',30,'MinPeakProminence',pProm);
            [valleysT locsT2]=findpeaks(-synchData.meanVLNAt,'MinPeakDistance',30,'MinPeakProminence',tProm);
            
            pTimes=synchData.time(locsP);
            tTimes=synchData.time(locsT);
            
            sPeakP=timetable(pTimes,peaksP);
            sPeakT=timetable(tTimes,peaksT);
            
            pS=withtol(pTimes,seconds(15));
            goodTimesT=sPeakT(pS,:);
            tS=withtol(tTimes,seconds(15));
            goodTimesP=sPeakP(tS,:);
            
            sortedOutP=ismember(tTimes,goodTimesT.tTimes);
            sortedOutT=ismember(pTimes,goodTimesP.pTimes);
            
            pTimes(find(sortedOutT==0))=[];
            tTimes(find(sortedOutP==0))=[];
            peaksP(find(sortedOutT==0))=[];
            peaksT(find(sortedOutP==0))=[];
            
            pTimes2=synchData.time(locsP2);
            tTimes2=synchData.time(locsT2);
            
            sValleysP=timetable(pTimes2,valleysP);
            sValleysT=timetable(tTimes2,valleysT);
            
            pS2=withtol(pTimes2,seconds(15));
            goodTimesT2=sValleysT(pS2,:);
            tS2=withtol(tTimes2,seconds(15));
            goodTimesP2=sValleysP(tS2,:);
            
            sortedOutP2=ismember(tTimes2,goodTimesT2.tTimes2);
            sortedOutT2=ismember(pTimes2,goodTimesP2.pTimes2);
            
            pTimes2(find(sortedOutT2==0))=[];
            tTimes2(find(sortedOutP2==0))=[];
            valleysP(find(sortedOutT2==0))=[];
            valleysT(find(sortedOutP2==0))=[];
            
            % Calculate shift in seconds
            etimePeaks=etime(datevec(pTimes),datevec(tTimes));
            etimeValleys=etime(datevec(pTimes2),datevec(tTimes2));
            
            etimeMean=round(mean(cat(1,etimePeaks,etimeValleys)));
            
            timeLagAll=cat(1,timeLagAll,etimePeaks,etimeValleys);
            
            timeNewSync=synchData.time+seconds(etimeMean);
            
            vlnaTimeNew=timetable(timeNewSync,synchData.meanVLNAt);
            
            synchData=sortrows(synchData,'time');
            synchData = rmmissing(synchData);
            vlnaTimeNew=sortrows(vlnaTimeNew,'timeNewSync');
            vlnaTimeNew = rmmissing(vlnaTimeNew);
            
            syncShift=synchronize(synchData,vlnaTimeNew);
            syncShift=rmmissing(syncShift);
            
            syncShift.Properties.VariableNames{'Var1_synchData'} = 'DBMVC';
            syncShift.Properties.VariableNames{'Var1_vlnaTimeNew'} = 'VLNAtemps';
            syncShift.meanVLNAt=[];
            plotSym='o';
            
            lnaTempsAll=cat(1,lnaTempsAll,syncShift.VLNAtemps);
        elseif inlist{ii,13}==2
            syncShift=synchData;
            syncShift.Properties.VariableNames{'Var1'} = 'DBMVC';
            syncShift.Properties.VariableNames{'meanVLNAt'} = 'VLNAtemps';
            plotSym='+';
        end
                     
        %% Plot dbmvc and vlna temp
        f3=figure('visible',dispFig);
        set(f3,'Position',[1500 500 1500 1000]);
        
        subplot(2,1,1)
        hold on;
        
        plot(time,dBmVCmean,'c');
        %plot(pTimes,peaksP,'o','MarkerfaceColor','b')
        %plot(pTimes2,-valleysP,'o','MarkerfaceColor','b')
        xlabel('Time [UTC]');
        ylabel('DBMVC [dB]');
        
        xlim([newTimes(1) newTimes(end)]);
        
        yyaxis right
        plot(timeTemp,VLnaTemp,'-b','linewidth',1.5)
        plot(timeTemp,meanVLNAt,'-r','linewidth',1.5)
        %plot(tTimes,peaksT,'o','MarkerfaceColor','b')
        %plot(tTimes2,-valleysT,'o','MarkerfaceColor','b')
        ylabel('Temperature [C]');
        
        formatOut = 'yyyymmdd_HHMMSS';
        timestring=datestr(time(1),formatOut);
        title([timestring],'interpreter','none');
        
        legend('DBMVC','VLNA Temperature','Mean VLNA Temp');
        xlim([newTimes(1) newTimes(end)]);
        
        subplot(2,1,2)
        hold on
        
        plot(syncShift.time,syncShift.DBMVC,'-c','linewidth',1.5);
        xlabel('Time [UTC]');
        ylabel('DBMVC [dB]');
        if ~isempty(syncShift.DBMVC)
            xlim([newTimes(1) newTimes(end)]);
        end
        
        yyaxis right
        
        plot(syncShift.time,syncShift.VLNAtemps,'-r','linewidth',1.5);
        ylabel('Temperature [C]');
        
        if ~isempty(syncShift.VLNAtemps)
            xlim([newTimes(1) newTimes(end)]);
        end
        legend('DBMVC sync','VLNA Temp sync');
        
        set(f3,'PaperPositionMode','auto')
        print(f3, [figdir 'lnaTemps_nsCal_' project '_',timestring,],'-dpng','-r0')
        
        %% Scatter plot vlna temp vs dbmvc
        f4=figure('visible',dispFig);
        set(f4,'Position',[1500 500 500 500]);
        
        hold on;
        
        scatter(syncShift.VLNAtemps,syncShift.DBMVC,plotSym,'MarkerFaceColor','b');
        
        xlabel('VLNA temperature [C]');
        ylabel('DBMVC [dB]');
        
        title([timestring],'interpreter','none');
        
        xlimits1=xlim;
        
        %fitAll1=polyfit(syncShift.VLNAtemps,syncShift.DBMVC,1);
        %fitAll1=orthDistRegression(syncShift.VLNAtemps,syncShift.DBMVC);
        fitOrth1=gmregress(syncShift.VLNAtemps,syncShift.DBMVC);
        fitAll1=[fitOrth1(2) fitOrth1(1)];
        xFit1 = xlimits1(1):0.1:xlimits1(2);
        yFit1 = polyval(fitAll1, xFit1);
        
        plot(xFit1, yFit1,'-r','linewidth',3);
                
        th=TextLocation(['y = ',num2str(fitAll1(1)),' x ',num2str(fitAll1(2))],'location','best');
        th.FontSize=11;
        th.Color='r';
        
        if inlist{ii,13}==1
            powTempSlope=cat(1,powTempSlope,fitAll1(1));
        end
        
        set(f4,'PaperPositionMode','auto')
        print(f4, [figdir 'lnaTempsScatter_nsCal_' project '_',timestring,],'-dpng','-r0')
        
        figure(f5);
        scatter(syncShift.VLNAtemps,syncShift.DBMVC,[],colorAll(ii,:),plotSym);
    end
    
    figure(f1);
    scatter(newTemp,powerDiffFile,[],colorAll(ii,:),plotSym);
    drawnow;
    
    figure(f6);
    scatter(newTemp,newDBMVC,[],colorAll(ii,:),plotSym);
    drawnow;
    
    legendCell{end+1}=timestring;
end
%% Scatter plot
legendCell{end+1}='Fit';

figure(f1)
xlimits=xlim;
ylimits=ylim;

fitOrth=gmregress(TempAll,powerDiff,1);
fitAll=[fitOrth(2) fitOrth(1)];
xFit = xlimits(1):0.1:xlimits(2);
yFit = polyval(fitAll, xFit);

plot(xFit, yFit,'-k','linewidth',3);

xlim(xlimits);
ylim(ylimits);

xlabel('Temperature [C]');
ylabel('Power difference [dB]');
title([project ' noise source calibration']);

textbp(['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
legend(legendCell,'interpreter','none','location','eastoutside');

set(f1,'PaperPositionMode','auto')
print(f1,[figdir 'nsCal_' project],'-dpng','-r0')  

if  any(inlist{:,13}==1) | any(inlist{ii,13}==2)
    figure(f5)
            
    xlabel('Temperature [C]');
    ylabel('DBMVC [dB]');
    title([project ' DBMVC vs VLNA temperature']);
    
    legend(legendCell,'interpreter','none','location','eastoutside');
    
    set(f5,'PaperPositionMode','auto')
    print(f5,[figdir 'nsCal_' project '_LNAtemps'],'-dpng','-r0')
end

figure(f6)
xlimitsD=xlim;
ylimitsD=ylim;

fitOrthD=gmregress(TempAll,dbmvcAll,1);
fitAllD=[fitOrthD(2) fitOrthD(1)];
xFitD = xlimitsD(1):0.1:xlimitsD(2);
yFitD = polyval(fitAllD, xFitD);

plot(xFitD, yFitD,'-k','linewidth',3);

xlim(xlimitsD);
ylim(ylimitsD);

xlabel('Temperature [C]');
ylabel('DBMVC [dB]');
title([project ' noise source calibration']);

textbp(['y = ',num2str(fitAllD(1)),' x + ',num2str(fitAllD(2))],'FontSize',14);
legend(legendCell,'interpreter','none','location','eastoutside');

set(f6,'PaperPositionMode','auto')
print(f6,[figdir 'nsCal_DBMVC_' project],'-dpng','-r0')

% Calculate means
if  any(inlist{:,13})==1
    disp(['Time lag between DBMVC and LNA temperature: ',num2str(mean(timeLagAll)),...
        ' +- ',num2str(std(timeLagAll)),' seconds']);
    disp(['Slope of DBMVC vs temperature data: ',num2str(mean(powTempSlope)),...
        ' +- ',num2str(std(powTempSlope))]);
    disp(['Mean LNA temperature for stable cases: ',num2str(mean(lnaTempsAll)),...
        ' +- ',num2str(std(lnaTempsAll)),' C']);
    
    meanTable=table([mean(timeLagAll);std(timeLagAll)],[mean(powTempSlope);std(powTempSlope)],...
        [nan;nan],[nan;nan],[nan;nan],[nan;nan],...
        'VariableNames',{'lnaTempLagSecs';' lnaGainChangePerC';'lnaTempReference';...
        'rxGainChangePerC';'podTempReference';'oceanScanBiasDB'});
    
  writetable(meanTable,[filedir,'tempsTable_',project,'.txt'],'Delimiter',' ');
end