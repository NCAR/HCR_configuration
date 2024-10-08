% Compare HCR temperatures

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc0';
freqData='10hz_combined';
qcVersion='';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'temperatures/'];


for jj=7:size(caseList,1)
    
    disp(['IOP ',num2str(jj),' of ',num2str(size(caseList,1))]);
    startTime=datetime(caseList(jj,1:6));
    endTime=datetime(caseList(jj,7:12));
    
    allFileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    eikT=[];
    polSwitchT=[];
    RfDetectorT=[];
    NST=[];
    timeTemp=[];
    VlnaT=[];
    HlnaT=[];
    XmitT=[];
    ploT=[];
    ps28VT=[];
    rotMotT=[];
    tiltMotT=[];
    insT=[];
    tailT=[];
    penFpgaT=[];
    penBoardT=[];
    rdsT=[];
    
    for ii=2:size(allFileList,2)
        tempFile=allFileList{ii};
        eikT=[eikT varFromCfRadialString(tempFile,'EikTemp')];
        polSwitchT=[polSwitchT varFromCfRadialString(tempFile,'PolarizationSwitchTemp')];
        RfDetectorT=[RfDetectorT varFromCfRadialString(tempFile,'RfDetectorTemp')];
        NST=[NST varFromCfRadialString(tempFile,'NoiseSourceTemp')];
        VlnaT=[VlnaT varFromCfRadialString(tempFile,'VLnaTemp')];
        HlnaT=[HlnaT varFromCfRadialString(tempFile,'HLnaTemp')];
        XmitT=[XmitT varFromCfRadialString(tempFile,'<XmitterTemp>')];
        ploT=[ploT varFromCfRadialString(tempFile,'PloTemp')];
        ps28VT=[ps28VT varFromCfRadialString(tempFile,'Ps28VTemp')];
        rdsT=[rdsT varFromCfRadialString(tempFile,'RdsInDuctTemp')];
        rotMotT=[rotMotT varFromCfRadialString(tempFile,'RotationMotorTemp')];
        tiltMotT=[tiltMotT varFromCfRadialString(tempFile,'TiltMotorTemp')];
        insT=[insT varFromCfRadialString(tempFile,'InsTemp')];
        tailT=[tailT varFromCfRadialString(tempFile,'TailconeTemp')];
        penFpgaT=[penFpgaT varFromCfRadialString(tempFile,'PentekFpgaTemp')];
        penBoardT=[penBoardT varFromCfRadialString(tempFile,'PentekBoardTemp')];
        
        startTimeIn=ncread(tempFile,'time_coverage_start')';
        startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
        timeTemp=[timeTemp startTimeFile];
    end
    
    %% Plot
    close all
    
    f1=figure('DefaultAxesFontSize',14);
    set(f1,'Position',[200 500 2000 1300]);
    set(f1,'renderer','painters');
    
    subplot(4,1,1)
    hold on;
    plot(timeTemp,penFpgaT,'linewidth',2);
    plot(timeTemp,penBoardT,'linewidth',2);
    plot(timeTemp,ps28VT,'linewidth',2);
    
    xlim([timeTemp(1) timeTemp(end)]);
    yticks(-100:10:100);
    ylim([10 80])
    ylabel('Temperature [C]');
    legend('PentekFpga','PentekBoard','Ps28V','Orientation','horizontal','location','northwest');
    grid on
    title([project,' IOP ',num2str(jj)]);
    
    subplot(4,1,2)
    hold on;
    plot(timeTemp,VlnaT,'linewidth',2);
    plot(timeTemp,HlnaT,'linewidth',2);
    plot(timeTemp,XmitT,'linewidth',2);
        
    xlim([timeTemp(1) timeTemp(end)]);
    yticks(-100:10:100);
    ylim([10 50])
    ylabel('Temperature [C]');
    legend('VLna','HLna','Xmitter','Orientation','horizontal','location','northwest');
    grid on
    
    subplot(4,1,3)
    hold on;
    plot(timeTemp,rotMotT,'linewidth',2);
    plot(timeTemp,tiltMotT,'linewidth',2);
    plot(timeTemp,eikT,'linewidth',2);
    plot(timeTemp,polSwitchT,'linewidth',2);
    plot(timeTemp,RfDetectorT,'linewidth',2);
    plot(timeTemp,NST,'linewidth',2);
    plot(timeTemp,ploT,'linewidth',2);
    plot(timeTemp,rdsT,'-k','linewidth',2);
    
    xlim([timeTemp(1) timeTemp(end)]);
    yticks(-100:10:100);
    ylim([0 60])
    ylabel('Temperature [C]');
    legend('RotationMotor','TiltMotor','Eik','PolarizationSwitch','RfDetector','NoiseSource','Plo','RdsInDuct','Orientation','horizontal','location','northwest');
    grid on
    
    subplot(4,1,4)
    hold on;
    plot(timeTemp,insT,'linewidth',2);
    plot(timeTemp,tailT,'linewidth',2);
    
    xlim([timeTemp(1) timeTemp(end)]);
    yticks(-100:10:100);
    ylim([0 50])
    ylabel('Temperature [C]');
    legend('Ins','Tailcone','Orientation','horizontal','location','northwest');
    grid on
        
    print([figdir,project,'_IOP',num2str(jj),'_',datestr(timeTemp(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timeTemp(end),'yyyymmdd_HHMMSS'),'_temperatures'],'-dpng','-r0')
    
end