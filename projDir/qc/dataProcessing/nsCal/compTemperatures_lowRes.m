% Compare HCR temperatures

clear all;
close all;

addpath('/h/eol/romatsch/git/private/utils/');

%figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/socrates/temperatures/'];
figdir=['/h/eol/romatsch/hcrCalib/nsCal/figs/qc1/otrec/temperatures/'];

lowResTempDir='/scr/snow2/rsfdata/projects/otrec/hcr/cfradial/moments/10hz/20190419/';
%lowResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/bad_ldr.10hz/20171113/'; % Lab cal 1
%lowResTempDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/labcal/10hz/20180921/'; % Lab cal 2

eikT=[];
polSwitchT=[];
RfDetectorT=[];
NST=[];
timeTemp=[];
VlnaT=[];
HlnaT=[];

allFileList=dir([lowResTempDir,'*.nc']);

for ii=1:size(allFileList,1)
    tempFile=[allFileList(ii).folder,'/',allFileList(ii).name];
    eikT=[eikT varFromCfRadialString(tempFile,'EikTemp')];
    polSwitchT=[polSwitchT varFromCfRadialString(tempFile,'PolarizationSwitchTemp')];
    RfDetectorT=[RfDetectorT varFromCfRadialString(tempFile,'RfDetectorTemp')];
    NST=[NST varFromCfRadialString(tempFile,'NoiseSourceTemp')];
    VlnaT=[VlnaT varFromCfRadialString(tempFile,'VLnaTemp')];
    HlnaT=[HlnaT varFromCfRadialString(tempFile,'HLnaTemp')];
    
    startTimeIn=ncread(tempFile,'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeTemp=[timeTemp startTimeFile];
end

f2=figure('DefaultAxesFontSize',14);
set(f2,'Position',[200 500 1500 500]);
hold on;
plot(timeTemp,eikT);
plot(timeTemp,polSwitchT);
plot(timeTemp,RfDetectorT);
plot(timeTemp,NST);
plot(timeTemp,VlnaT);
plot(timeTemp,HlnaT);
legend('EikTemp','PolSwitchTemp','RfDetTemp','NoisSourceTemp','VLnaTemp','HLnaTemp','location','best');

ylabel('Temperature [C]');

formatOut = 'yyyymmdd_HHMMSS';
timestring=datestr(timeTemp(1),formatOut);
timestring2=datestr(timeTemp(end),formatOut);
title([timestring,' to ',timestring2],'interpreter','none');

set(gcf,'PaperPositionMode','auto')
print(f2, [figdir,'temperatures_',timestring,'_to_',timestring2],'-dpng','-r0');

disp(['Mean VLNA temperature is ',num2str(nanmean(VlnaT)),' C']);