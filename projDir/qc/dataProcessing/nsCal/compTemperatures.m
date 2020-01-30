% Compare HCR temperatures

clear all;
close all;

project='socrates';

addpath('/h/eol/romatsch/gitPriv/utils/');

if strcmp(project,'socrates')
    figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc2/socrates/temperatures/'];
    highResTempDir='/scr/snow2/rsfdata/projects/socrates/hcr/qc/temperatures1s/';
elseif strcmp(project,'cset')
    figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/cset/'];
    highResTempDir='/h/eol/romatsch/data/hcrCalib/temps/';
end

allFileList=dir([highResTempDir,'*.txt']);

for jj=1:size(allFileList,1)
    close all
    if strcmp(project,'socrates')
        tempFile=[allFileList(jj).folder,'/',allFileList(jj).name];
        indata=txtTable2matTable(tempFile,',');
    elseif strcmp(project,'cset')
        tempFile=[highResTempDir,'cset_temps.txt'];
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
    
    f2=figure('DefaultAxesFontSize',14);
    set(f2,'Position',[200 500 1500 500]);
    
    hold on;
    plot(timeTemp,EikTemp);
    plot(timeTemp,PolSwitchTemp);
    plot(timeTemp,RfDetTemp);
    plot(timeTemp,NoisSourceTemp);
    plot(timeTemp,VLnaTemp);
    plot(timeTemp,HLnaTemp);
    
    ylabel('Temperature [C]');
    
    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(timeTemp(1),formatOut);
    timestring2=datestr(timeTemp(end),formatOut);
    title([timestring,' to ',timestring2],'interpreter','none');
    
    legend('EikTemp','PolSwitchTemp','RfDetTemp','NoisSourceTemp','VLnaTemp','HLnaTemp','location','best');
    
    set(gcf,'PaperPositionMode','auto')
    print(f2, [figdir,'temperatures_',timestring,'_to_',timestring2],'-dpng','-r0')
end