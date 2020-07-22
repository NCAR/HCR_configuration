function [refTemp refTempStd] = f_getRefTemp(indir,startTime,endTime,whichTemp)
%Get reference pod temperature from the lab bench top calibration
allFiles=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',0);

eikT=[];
polSwitchT=[];
RfDetectorT=[];
NST=[];
vlnaT=[];
hlnaT=[];
timeTemp=[];

for ii=1:size(allFiles,2)
    eikT=[eikT varFromCfRadialString(allFiles{ii},'EikTemp')];
    polSwitchT=[polSwitchT varFromCfRadialString(allFiles{ii},'PolarizationSwitchTemp')];
    RfDetectorT=[RfDetectorT varFromCfRadialString(allFiles{ii},'RfDetectorTemp')];
    NST=[NST varFromCfRadialString(allFiles{ii},'NoiseSourceTemp')];
    vlnaT=[vlnaT varFromCfRadialString(allFiles{ii},'VLnaTemp')];
    hlnaT=[hlnaT varFromCfRadialString(allFiles{ii},'HLnaTemp')];
      
    startTimeIn=ncread(allFiles{ii},'time_coverage_start')';
    startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
        str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
    timeTemp=[timeTemp startTimeFile];
end

eikT(eikT>500)=nan;
polSwitchT(polSwitchT>500)=nan;
RfDetectorT(RfDetectorT>500)=nan;
NST(NST>500)=nan;
vlnaT(vlnaT>500)=nan;
hlnaT(hlnaT>500)=nan;

meanTemp=nanmean([eikT;polSwitchT;RfDetectorT;NST],1);

figure;
set(gcf,'Position',[200 500 1200 600]);
hold on;
plot(timeTemp,meanTemp);
plot(timeTemp,eikT);
plot(timeTemp,polSwitchT);
plot(timeTemp,RfDetectorT);
plot(timeTemp,NST);
plot(timeTemp,vlnaT);
plot(timeTemp,hlnaT);
legend('meanTemp','eikTemp','polSwitchTemp','RfDetectorTemp','NoiseSourceTemp','VlnaTemp','HlnaTemp');

if strcmp(whichTemp,'podTemp')
    refTemp=mean(meanTemp);
    refTempStd=std(meanTemp);
elseif strcmp(whichTemp,'vlnaTemp')
    refTemp=mean(vlnaT);
    refTempStd=std(vlnaT);
else
    disp('Requested temperature must be podTemp or vlnaTemp.');
end
end

