% When noise cal was done with high resolution temperature data, run this
% code to refine the results

clear all;
close all;

project='cset';

addpath('/h/eol/romatsch/gitPriv/utils/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/NSCAL/functions/');

figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc2/' project '/podTemps/'];

if strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; %socrates
    highResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
elseif strcmp(project,'cset')
    indir='/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/'; %cset raw
    highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
elseif strcmp(project,'aristo')
    indir='/scr/snow2/rsfdata/projects/aristo-17/hcr/cfradial/moments/10hz/'; %aristo raw
    highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
elseif strcmp(project,'otrec')
    indir='/scr/snow1/rsfdata/projects/otrec/hcr/qcfield/cfradial/moments/10hz/'; % otrec field
    highResTempDir='/scr/snow1/rsfdata/projects/otrec/hcr/qc0/temperatures/';
else
    disp('Project name not valid.')
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

meanTable=readtable([filedir,'tempsTable_',project,'.txt']);

slope=meanTable.lnaGainChangePerC(1);
%slope=0.175;

if strcmp(project,'socrates')
    [refTemp refTempStd]=f_getRefTemp([indir '20171113/'],datetime(2017,11,13,23,11,30),datetime(2017,11,13,23,32,17),'vlnaTemp');
elseif strcmp(project,'cset') | strcmp(project,'aristo')
    %[refTemp refTempStd]=f_getRefTemp(['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/20140825/'],datetime(2014,8,25,22,30,0),datetime(2014,8,25,23,0,0),'vlnaTemp');
    [refTemp refTempStd]=f_getRefTemp(['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/20160822/'],datetime(2016,8,22,21,0,0),datetime(2016,8,22,22,0,0),'vlnaTemp');
elseif strcmp(project,'otrec')
    [refTemp refTempStd]=f_getRefTemp(['/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/20190619/'],datetime(2019,6,19,20,0,0),datetime(2019,6,19,22,0,0),'vlnaTemp');
end

%% Loop through cases

powerDiff=[];
TempAll=[];
dbmvcAll=[];

f1=figure;
set(f1,'Position',[200 100 1200 1000]);
hold on

f2=nan;
f3=nan;
f4=nan;

f5=figure;
set(f5,'Position',[200 100 1200 1000]);
hold on
legendCell={};

timeLagAll=[];
powTempSlope=[];
lnaTempsAll=[];

f6=figure;
set(f6,'Position',[200 100 1200 1000]);
hold on
legendCell={};

colorAll=jet(size(inlist,1));

for ii=1:size(inlist,1)
    disp(['Case ' num2str(ii) ' from ' num2str(size(inlist,1))]);
    startTime=datetime(inlist{ii,1:6});
    endTime=datetime(inlist{ii,7:12});
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if inlist{ii,13}==3
        continue
    end
    
    if inlist{ii,13}==1
        plotSym='o';
    else
        plotSym='+';
    end
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
    
    if strcmp(project,'socrates')
        tempFile=highResTempFiles_socrates(startTime,endTime,highResTempDir);
        indata=txtTable2matTable(tempFile,',');
    elseif strcmp(project,'cset') | strcmp(project,'aristo') | strcmp(project,'otrec')
        if ii==1
            tempFile=[highResTempDir,project,'_temps.txt'];
            tempnames={'count','year','month','day','hour','min','sec','unix_time',...
                'unix_day','XmitterTemp','PloTemp','EikTemp','VLnaTemp','HLnaTemp',...
                'PolarizationSwitchTemp','RfDetectorTemp','NoiseSourceTemp','Ps28VTemp',...
                'RdsInDuctTemp','RotationMotorTemp','TiltMotorTemp','CmigitsTemp',...
                'TailconeTemp','PentekFpgaTemp','PentekBoardTemp'};
            indata=readtable(tempFile);
            indata.Properties.VariableNames=tempnames;
        end
    end
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Smooth LNA and pod temperature data and adjust for time shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meanVLNAt=movmean(VLnaTemp,20);
    meanPodT=movmean(Temp,80);
    powersMean=movmean(dBmVCmean',20);
    powers=timetable(time,powersMean);
    vlnatemps=timetable(timeTemp,meanVLNAt,meanPodT);
    
    synchData=synchronize(powers,vlnatemps,'secondly','mean');
    
    synchData=rmmissing(synchData);
    
    timeNewSync=synchData.time+seconds(round(meanTable.lnaTempLagSecs(1)));
    
    vlnaTimeNew=timetable(timeNewSync,synchData.meanVLNAt);
    
    syncShift=synchronize(synchData,vlnaTimeNew);
    syncShift=rmmissing(syncShift);
    
    syncShift.meanVLNAt=[];
    syncShift.Properties.VariableNames{'Var1'} = 'meanVLNAt';
    
    % Adjust power data with LNA temperature data
    degDiff=refTemp-syncShift.meanVLNAt;
    gainChange=slope.*degDiff;
    
    newDBMVC=syncShift.powersMean+gainChange;
    
    %Calculate theoretical DBMVC value
    Bn=1/pulseWidth(1)*0.7;
    
    ENR_Corr=10^(ENR_Quinstar/10) + (T0-(syncShift.meanPodT + 273.15))./T0;
    
    Ts_ENR_Corr = ENR_Corr.*T0+T0;
    Pn_ENR_Corr = K.*Ts_ENR_Corr.*Bn;
    
    Pn_ENR_Corr_dB=10.*log10(Pn_ENR_Corr);
    
    powerDiffFile=newDBMVC-(Pn_ENR_Corr_dB+30);
    
    if  ishghandle(f2)
        close(f2);
    end
    
    f2=figure;
    set(f2,'Position',[1500 500 1500 500]);
    
    hold on;
    
    plot(time,dBmVCmean,'c');
    plot(syncShift.time,newDBMVC,'-b','linewidth',2);
    plot(syncShift.time,Pn_ENR_Corr_dB+30,'-r','linewidth',2);
    xlabel('Time [UTC]');
    ylabel('DBMVC [dB]');
    xlim([time(1) time(end)]);
    
    yyaxis right
    plot(timeTemp,Temp,'-g','MarkerFaceColor','g')
    plot(syncShift.time,syncShift.meanPodT,'-k','linewidth',2)
    ylabel('Temperature [C]');
    
    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(time(1),formatOut);
    title([timestring],'interpreter','none');
    
    legend('DBMVC raw','DBMVC corr','Noise Source','Pod temp raw','Pod temp smooth');
    
    
    set(f2,'PaperPositionMode','auto')
    print(f2, [figdir 'nsCal_' project '_',timestring],'-dpng','-r0')
        
    figure(f5);
    scatter(syncShift.meanVLNAt,syncShift.powersMean,[],colorAll(ii,:),plotSym);
    
    figure(f1);
    scatter(syncShift.meanPodT,powerDiffFile,[],colorAll(ii,:),plotSym);
    drawnow;
    
    figure(f6);
    scatter(syncShift.meanPodT,newDBMVC,[],colorAll(ii,:),plotSym);
    drawnow;
    
    legendCell{end+1}=timestring;
    
    powerDiff=cat(1,powerDiff,powerDiffFile);
    TempAll=cat(1,TempAll,syncShift.meanPodT);
    dbmvcAll=cat(1,dbmvcAll,newDBMVC);
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
ylabel('DBMVC minus noise source [dB]');
title([project ' noise source calibration']);

textbp(['y = ',num2str(fitAll(1)),' x + ',num2str(fitAll(2))],'FontSize',14);
legend(legendCell,'interpreter','none','location','eastoutside');

set(f1,'PaperPositionMode','auto')
print(f1,[figdir 'nsCal_' project],'-dpng','-r0')

figure(f5)

xlabel('Temperature [C]');
ylabel('DBMVC [dB]');
title([project ' DBMVC vs VLNA temperature']);

legend(legendCell,'interpreter','none','location','eastoutside');

set(f5,'PaperPositionMode','auto')
print(f5,[figdir 'nsCal_' project '_LNAtemps'],'-dpng','-r0')

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

if strcmp(project,'socrates')
    [podTemp podTempStd]=f_getRefTemp([indir '20171113/'],datetime(2017,11,13,23,11,30),datetime(2017,11,13,23,32,17),'podTemp');
elseif strcmp(project,'cset') | strcmp(project,'aristo')
    %[podTemp podTempStd]=f_getRefTemp(['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/20140825/'],datetime(2014,8,25,22,30,0),datetime(2014,8,25,23,0,0),'podTemp');
    [podTemp podTempStd]=f_getRefTemp(['/scr/snow2/rsfdata/projects/cset/hcr/cfradial/moments/10hz/20160822/'],datetime(2016,8,22,21,0,0),datetime(2016,8,22,22,0,0),'podTemp');
elseif strcmp(project,'otrec')
    [podTemp podTempStd]=f_getRefTemp(['/scr/snow1/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/20190619/'],datetime(2019,6,19,20,0,0),datetime(2019,6,19,22,0,0),'podTemp');
end

meanTable.lnaTempReference(1)=refTemp;
meanTable.lnaTempReference(2)=refTempStd;
meanTable.rxGainChangePerC(1)=fitAll(1);
meanTable.podTempReference(1)=podTemp;
meanTable.podTempReference(2)=podTempStd;

writetable(meanTable,[filedir,'tempsTable_',project,'.txt'],'Delimiter',' ');