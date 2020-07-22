% Compare HCR temperatures

clear all;
close all;

project='cset';

addpath('/h/eol/romatsch/gitPriv/utils/');

if strcmp(project,'socrates')
    figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc2/socrates/temperatures/'];
    highResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/data/socrates/temperatures/';
    allFileList=dir([highResTempDir,'*.txt']);
elseif strcmp(project,'cset') | strcmp(project,'aristo')
    figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc2/',project,'/temperatures/'];
    if strcmp(project,'cset')
        highResTempDir='/scr/snow2/rsfdata/projects/cset/hcr/txt/';
    elseif strcmp(project,'aristo')
        highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
    end
    allFileList=1;
end

for jj=1:size(allFileList,1)
    close all
    if strcmp(project,'socrates')
        tempFile=[allFileList(jj).folder,'/',allFileList(jj).name];
        indata=txtTable2matTable(tempFile,',');
    elseif strcmp(project,'cset') | strcmp(project,'aristo')
        if strcmp(project,'cset')
            tempFile=[highResTempDir,'CSET.temperatures.txt'];
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
        
    EikTemp=indata.EikTemp;
    PolSwitchTemp=indata.PolarizationSwitchTemp;
    RfDetTemp=indata.RfDetectorTemp;
    NoisSourceTemp=indata.NoiseSourceTemp;
    VLnaTemp=indata.VLnaTemp;
    HLnaTemp=indata.HLnaTemp;
        
    timeTemp=datetime(indata.year,indata.month,indata.day,indata.hour,indata.min,indata.sec);
    
    % Find gaps between flights
    threshold = hours(3);
    gap = diff(timeTemp);
    idx = find(gap > threshold);
    idx=[1;idx;length(timeTemp)];
    
    for ii=1:length(idx)-1
        close all
        f2=figure('DefaultAxesFontSize',14);
        set(f2,'Position',[200 500 1500 500]);
        
        hold on;
        plot(timeTemp(idx(ii)+1:idx(ii+1)),EikTemp(idx(ii)+1:idx(ii+1)));
        plot(timeTemp(idx(ii)+1:idx(ii+1)),PolSwitchTemp(idx(ii)+1:idx(ii+1)));
        plot(timeTemp(idx(ii)+1:idx(ii+1)),RfDetTemp(idx(ii)+1:idx(ii+1)));
        plot(timeTemp(idx(ii)+1:idx(ii+1)),NoisSourceTemp(idx(ii)+1:idx(ii+1)));
        plot(timeTemp(idx(ii)+1:idx(ii+1)),VLnaTemp(idx(ii)+1:idx(ii+1)));
        plot(timeTemp(idx(ii)+1:idx(ii+1)),HLnaTemp(idx(ii)+1:idx(ii+1)));
        
        ylabel('Temperature [C]');
        
        formatOut = 'yyyymmdd_HHMMSS';
        timestring=datestr(timeTemp(idx(ii)+1),formatOut);
        timestring2=datestr(timeTemp(idx(ii+1)),formatOut);
        title([timestring,' to ',timestring2],'interpreter','none');
        
        legend('EikTemp','PolSwitchTemp','RfDetTemp','NoisSourceTemp','VLnaTemp','HLnaTemp','location','best');
        
        set(gcf,'PaperPositionMode','auto')
        print(f2, [figdir,'temperatures_',timestring,'_to_',timestring2],'-dpng','-r0')
    end
end