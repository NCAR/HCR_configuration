% Calculate PID from HCR HSRL combined data

clear all
%close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='2hzMerged'; % 10hz, 100hz, or 2hz

startTime=datetime(2018,1,19,4,20,0); %Wang_Rauber
endTime=datetime(2018,1,19,4,40,0); %Wang_Rauber

ylimits=[0 4];

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/pid_hcr/',project,'/'];

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if ~isempty(fileList)
    %% Load data
    
    %HCR data
    data.PID_COMBINED=[];
    data.PID_HCR=[];
    data.PID_HSRL=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    data.asl=data.asl./1000;
    
    %% Scales and units
    cscale_hsrl=[0,0,1.0;0,1,0;1,0.67,0;1,0,1;0,1,1;1,0.67,0];
    cscale_hcr=[0,0,1.0; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
    cscale_comb=[0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
    
    units_str_hsrl={'Cloud liquid','Drizzle',...
        'Aerosol1','SLW','Ice crystals','Aerosol2'};
    units_str_hcr={'Cloud liquid','Drizzle',...
        'Rain','SLW','Ice crystals','Snow','Wet snow/rimed ice'};
    units_str_comb={'Cloud liquid','Drizzle','Rain',...
        'SLW','Ice crystals','Snow','Wet snow/rimed ice','Aerosols'};
    
    %% Plot
    close all
    f1=figure('DefaultAxesFontSize',12,'Position',[400 300 1200 900]);
    
    s1=subplot(3,1,1);
    surf(data.time,data.asl,data.PID_HSRL,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([.5 6.5]);
    colormap(s1,cscale_hsrl);
    cb=colorbar;
    cb.Ticks=1:6;
    cb.TickLabels=units_str_hsrl;
    ylabel('Altitude (km)');
    title(['HSRL particle ID']);
    
    s2=subplot(3,1,2);
    surf(data.time,data.asl,data.PID_HCR,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([.5 7.5]);
    colormap(s2,cscale_hcr);
    cb=colorbar;
    cb.Ticks=1:7;
    cb.TickLabels=units_str_hcr;
    ylabel('Altitude (km)');
    title(['HCR particle ID']);
    
    s3=subplot(3,1,3);
    surf(data.time,data.asl,data.PID_COMBINED,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([.5 8.5]);
    colormap(s3,cscale_comb);
    cb=colorbar;
    cb.Ticks=1:8;
    cb.TickLabels=units_str_comb;
    ylabel('Altitude (km)');
    title(['Combined particle ID']);
    
    %     f4=figure('DefaultAxesFontSize',12,'Position',[400 300 1300 500]);
    %
    %     s1=subplot(1,1,1);
    %     fig1=surf(data.time,data.asl,data.PID+1,'edgecolor','none');
    %     view(2);
    %     ylim(ylimits);
    %     %xlim([data.time(1),data.time(end)]);
    %     xlim([datetime(2018,1,29,1,50,0),data.time(end)]);
    %     caxis([.5 9.5]);
    %     colormap(s1,cscale_comb);
    %     cb=colorbar;
    %     cb.Ticks=1:9;
    %     cb.TickLabels=units_str_comb;
    %     ylabel('Altitude (km)');
    %     title(['Particle ID Combined']);
    
end