% Perform noise cal

clear all;
close all;

project='socrates';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/h/eol/romatsch/papers/HCRcalibration/figs/'];

indir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; %socrates
highResTempDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/temperatures1s/';

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

for ii=24:24
    
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
    
    %Get high resolution temperature data if wanted
    if inlist{ii,13}==1 | inlist{ii,13}==2
        if strcmp(project,'socrates')
            tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);
            indata=txtTable2matTable(tempFile,',');
        elseif strcmp(project,'cset') | strcmp(project,'aristo') | strcmp(project,'otrec')
            if ii==1
                if strcmp(project,'cset')
                    tempFile=[highResTempDir,'CSET.temperatures.txt'];
                elseif strcmp(project,'aristo') | strcmp(project,'otrec')
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
    nanInd=find(isnan(dBmVCmean) | isnan(dnTime));
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % When using high resolution temperature data, plot LNA temperature vs
    % powers to get that correlation too
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Smooth LNA temperature data
    meanVLNAt=movmean(VLnaTemp,20);
    powers=timetable(time,dBmVCmean');
    vlnatemps=timetable(timeTemp,meanVLNAt);
    
    synchData=synchronize(powers,vlnatemps,'secondly','mean');
    
    synchData=rmmissing(synchData);
    
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
    
    % Adjust power data with LNA temperature data
    fitOrth1=gmregress(syncShift.VLNAtemps,syncShift.DBMVC);
    fitAll1=[fitOrth1(2) fitOrth1(1)];
    
    degDiff=mean(syncShift.VLNAtemps)-syncShift.VLNAtemps;
    gainChange=fitOrth1(2).*degDiff;
    newDBMVC=syncShift.DBMVC+gainChange;
    
    %% Plot dbmvc and vlna temp
    close all
    
    wi=10;
    hi=4;
    
    fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi],'renderer','painters');
    fig1.PaperPositionMode = 'manual';
    fig1.PaperUnits = 'inches';
    fig1.Units = 'inches';
    fig1.PaperPosition = [0, 0, wi, hi];
    fig1.PaperSize = [wi, hi];
    fig1.Resize = 'off';
    fig1.InvertHardcopy = 'off';
    
    set(fig1,'color','w');
    
    %%%%%%%%%%%%%%%%%%%%%%%% VEL_RAW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,3,1:2);
    hold on;
    outerpos1 = ax1.Position;
    ax1.Position = [outerpos1(1)-0.07 outerpos1(2)+0.02 outerpos1(3)+0.02 outerpos1(4)+0.02];
    
    plot(time,dBmVCmean,'r');
    ylabel('DBMVC (dB)');
    
    xlim([newTimes(1) newTimes(end)]);
    ylim([-89.65 -89.05]);
    text(datetime(2018,1,15,22,2,5),-89.05,'(a)','fontweight','bold','fontsize',13);
    
    yyaxis right
    plot(timeTemp,VLnaTemp,'-c','linewidth',1.)
    plot(timeTemp,meanVLNAt,'-b','linewidth',1.5)
    ylabel('VLNA temperature (C)');
    ax = gca;
    ax.YColor = 'b';
    
    leg=legend('DBMVC averaged','VLNA temp.','VLNA temp. smoothed',...
        'orientation','horizontal','fontsize',10,'Position',[0.13,0.91,0.4098,0.0547]);
    leg.ItemTokenSize = [15,18];
    xlim([newTimes(1) newTimes(end)]);
    ylim([24.5 30.5]);
    
    grid on
    
    ax2=subplot(2,3,4:5);
    hold on;
    outerpos1 = ax2.Position;
    ax2.Position = [outerpos1(1)-0.07 outerpos1(2)+0.0 outerpos1(3)+0.02 outerpos1(4)+0.02];
    
    plot(syncShift.time,syncShift.DBMVC,'-r','linewidth',1.5);
    plot(syncShift.time,newDBMVC,'-m','linewidth',1.5);
    ylabel('DBMVC (dB)');
    xlim([newTimes(1) newTimes(end)]);
    ylim([-89.65 -89.05]);
    text(datetime(2018,1,15,22,2,5),-89.05,'(b)','fontweight','bold','fontsize',13);
    
    yyaxis right
    grid on
    
    plot(syncShift.time,syncShift.VLNAtemps,'-b','linewidth',1.5);
    ylabel('VLNA temperature (C)');
    ax = gca;
    ax.YColor = 'b';
    
    xlim([newTimes(1) newTimes(end)]);
    ylim([24.5 30.5]);
    
    leg=legend('DBMVC resampled','DBMVC corr.','VLNA temp. shifted',...
        'orientation','horizontal','fontsize',10,'Position',[0.0951,0.415,0.4677,0.0547]);
    leg.ItemTokenSize = [15,18];
    
    ax3=subplot(2,3,3);
    hold on;
    outerpos1 = ax3.Position;
    ax3.Position = [outerpos1(1)+0.0 outerpos1(2)-0.44 outerpos1(3)+0.07 outerpos1(4)+0.47];
    
    scatter(syncShift.VLNAtemps,syncShift.DBMVC,plotSym,'MarkerFaceColor','b');
    
    xlabel('VLNA temperature (C)');
    ylabel('DBMVC (dB)');
    
    xlimits1=xlim;
    text(26.1,-89.13,'(c)','fontweight','bold','fontsize',13);
        
    xFit1 = xlimits1(1):0.1:xlimits1(2);
    yFit1 = polyval(fitAll1, xFit1);
    
    plot(xFit1, yFit1,'-r','linewidth',3);
    
    print(fig1, [figdir 'lnaTemps.png'],'-dpng','-r0')
    
end
