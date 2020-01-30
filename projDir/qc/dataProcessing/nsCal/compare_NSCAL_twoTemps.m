% When noise cal was done with high resolution temperature data, run this
% code to refine the results

clear all;
close all;

project='socrates';

addpath('/h/eol/romatsch/git/private/utils/');
addpath('/h/eol/romatsch/git/private/process_HCR/NSCAL/functions/');

figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/' project '/refined/compare/'];

if strcmp(project,'socrates')
    indir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; %socrates
    highResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/data/socrates/temperatures1s/';
    newTempFile='/scr/rain1/rsfdata/projects/socrates/hcr/qc/data/socrates/tempVsGain/deltaGains.txt';
else
    disp('Project name must be socrates.')
    return
end

filedir='/h/eol/romatsch/hcrCalib/nsCal/inFiles/';
infile=['cal_' project '.dat'];

inlist=readtable([filedir infile]);
%Infile has start and end dates. In last column it has 0, 1, or 2.
% 0: Use low res temperature data directly from cfradial file
% 1: Use high resolution temperature data when LNA heater is working well,
% i.e. it cycles in regular waves. LNA temperature data will be adjusted for
% time delay between LNA temperatues and power waves.
% 2: LNA heaters are not working properly. High resolution temperature data
% will be used but LNA temperatures will not be adjusted for time lag.

T0=290;
K=1.38e-23;
ENR_Quinstar = 20.84; %dB

%get new temperature data
fileID = fopen(newTempFile,'r');
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
adjTemp = textscan(fileID,formatSpec,'HeaderLines',25);
fclose(fileID);

adjTempTime=datetime(adjTemp{1,2:7});
adjLNAt=adjTemp{1,12};

%% Loop through cases
%
% powerDiff=[];
% TempAll=[];
% dbmvcAll=[];

f3=nan;

if any(inlist{:,13})==1 | any(inlist{:,13})==2
   
    timeLagAll=[];
    powTempSlope=[];
    lnaTempsAll=[];
end

colorAll=jet(size(inlist,1));
%plotSym='^';

for ii=24:size(inlist,1)
    disp(['Case ' num2str(ii) ' from ' num2str(size(inlist,1))]);
    startTime=datetime(inlist{ii,1:6});
    endTime=datetime(inlist{ii,7:12});
    
    fileList=makeFileList(indir,startTime,endTime,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get all relevant data and remove times that are not wanted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dBmVC=[];
    range=[];
    time=[];
    Temp=[];
    timeTemp=[];
    pulseWidth=[];
    
    for jj=1:size(fileList,2)
        disp(['File ',num2str(jj),' from ',num2str(size(fileList,2))]);
        dBmVC=cat(2,dBmVC,ncread(fileList{jj},'DBMVC'));
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
    time(tooSmall)=[];
    
    % Get rid of data that is not within the wanted time interval
    outOfBoundsInds=find(time<startTime | time>endTime);
    
    dBmVC(:,outOfBoundsInds)=[];
    dBmVCmean=nanmean(dBmVC,1);
    
    time(outOfBoundsInds)=[];
        
    %Get high resolution temperature data if wanted
    if inlist{ii,13}==1 | inlist{ii,13}==2
       if strcmp(project,'socrates')
           tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);
       end
       indata=txtTable2matTable(tempFile,',');
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
       
    %Interpolate Temperature data
    nanIndT=find(isnan(Temp) | isnan(dnTimeTemp));
    Temp(nanIndT)=[];
    dnTimeTemp(nanIndT)=[];
    timeTemp(nanIndT)=[];
    fitTemp=polyfit(dnTimeTemp,Temp,1);    
    newTemp=polyval(fitTemp,DNnewTimes);
    
    if  ishghandle(f3)
        close(f3);
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When using high resolution temperature data, plot LNA temperature vs
    % powers to get that correlation too
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                     
        % Plot dbmvc and vlna temp
        f3=figure;
        set(f3,'Position',[200 100 1500 1000]);
        
        subplot(2,1,1)
        hold on;
        
        plot(time,dBmVCmean,'c');
        xlabel('Time [UTC]');
        ylabel('DBMVC [dB]');
        
        xlim([newTimes(1) newTimes(end)]);
        
        yyaxis right
        plot(timeTemp,VLnaTemp,'-b','linewidth',1.5)
        plot(timeTemp,meanVLNAt,'-r','linewidth',1.5)
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
        xlim([newTimes(1) newTimes(end)]);
        
        yyaxis right
   
        plot(syncShift.time,syncShift.VLNAtemps,'-r','linewidth',1.5);
        ylabel('Temperature [C]');
        
        xlim([newTimes(1) newTimes(end)]);
        
        % get new temperature data at right time
        adjInd=find(adjTempTime>=newTimes(1) & adjTempTime<=newTimes(end));
        plot(adjTempTime(adjInd),adjLNAt(adjInd),'-g','linewidth',1.5);
       
       
         legend('DBMVC sync','VLNA Temp sync','VLNA Temp database');
        
        set(f3,'PaperPositionMode','auto')
        print(f3, [figdir 'compare_lnaTemps_nsCal_' project '_',timestring,],'-dpng','-r0')
               
    end
end
