% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

b_drizz = 0.52; % Z<-17 dBZ
b_rain = 0.68; % Z>-17 dBZ
%alpha = 0.21;

ylimUpper=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/scr/sci/romatsch/liquidWaterHCR/'];
figdir=['/home/romatsch/plots/HCR/liquidWater/clearVScloudy/'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/liquidWater_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    %data.U_SURF=[];
    %data.V_SURF=[];
    %data.SST=[];
    %data.TEMP=[];
    %data.PRESS=[];
    %data.RH=[];
    data.TOPO=[];
    data.FLAG=[];   
    data.DBMVC=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    frq=ncread(fileList{1},'frequency');
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    %% Clear vs cloudy
    
    % Find ocean surface gate
    [linInd maxGate rangeToSurf] = hcrSurfInds(data);
    
    % Reflectivity
    reflTemp=data.DBZ;
    
    % Remove bang
    reflTemp(data.FLAG==6)=nan;
    
    reflLin=10.^(reflTemp./10);
    reflOceanLin=nan(size(data.time));
    reflNoOceanLin=nan(size(data.time));
    
    for ii=1:length(data.time)
        if (~(maxGate(ii)<10 | maxGate(ii)>size(reflLin,1)-5)) & ~isnan(maxGate(ii))
            reflRay=reflLin(:,ii);
            reflOceanLin(ii)=sum(reflRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
            reflNoOceanLin(ii)=sum(reflRay(1:maxGate(ii)-6),'omitnan');
        end
    end
    
    reflFrac=reflNoOceanLin./reflOceanLin;
    
    % Power
    powerTemp=data.DBMVC;
    
    % Remove bang
    powerTemp(data.FLAG==6)=nan;
    
    powerLin=10.^(powerTemp./10);
    powerOceanLin=nan(size(data.time));
    powerNoOceanLin=nan(size(data.time));
    
    for ii=1:length(data.time)
        if (~(maxGate(ii)<10 | maxGate(ii)>size(powerLin,1)-5)) & ~isnan(maxGate(ii))
            powerRay=powerLin(:,ii);
            powerOceanLin(ii)=sum(powerRay(maxGate(ii)-5:maxGate(ii)+5),'omitnan');
            powerNoOceanLin(ii)=sum(powerRay(1:maxGate(ii)-6),'omitnan');
        end
    end
    
    powerFrac=powerNoOceanLin./powerOceanLin;
    
    %% Plot
    close all
        
    f1 = figure('Position',[200 500 1500 1100],'DefaultAxesFontSize',12);
    
    s1=subplot(5,1,3);
    hold on
    l0=plot(data.time,reflFrac,'-r','linewidth',1);
    l1=plot(data.time,powerFrac,'-b','linewidth',1);
    ylabel('Fraction');
    ylim([0 0.001]);
    grid on
    
    xlim([data.time(1),data.time(end)]);
    
    legend([l0 l1],{'Reflectivity','Power'},'location','northeast');
    title({'Atmosphere divided by ocean, power or reflectivity'})
    s1pos=s1.Position;
    
    s3=subplot(5,1,4);
    l2=plot(data.time,reflFrac,'-r','linewidth',1);
    ylabel('Fraction');
    ylim([0 0.00001]);
    grid on
    set(gca,'YColor','k');
    
    xlim([data.time(1),data.time(end)]);
    
    legend({'Reflectivity'},'location','northeast');
    title('Atmosphere divided by ocean, reflectivity zoomed in');
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    s4=subplot(5,1,2);
    hold on
    l0=plot(data.time,reflOceanLin,'-r','linewidth',1);
    ylabel('Reflectivity');
    %ylim([0 200000]);
    
    yyaxis right
    l1=plot(data.time,powerOceanLin,'-b','linewidth',1);
    ylabel('Power');
    %ylim([0 0.0005]);
    grid on
    set(gca,'YColor','k');
    
    xlim([data.time(1),data.time(end)]);
    
    legend([l0 l1],{'Reflectivity','Power'},'location','northeast');
    title('Linear ocean surface power and reflectivity');
    
    s4pos=s4.Position;
    s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
    
    s5=subplot(5,1,1);
    hold on
    l0=plot(data.time,reflNoOceanLin,'-r','linewidth',1);
    ylabel('Reflectivity');
    ylim([0 0.2]);
    
    yyaxis right
    l1=plot(data.time,powerNoOceanLin,'-b','linewidth',1);
    ylabel('Power');
    ylim([0 0.00000002]);
    grid on
    set(gca,'YColor','k');
    
    xlim([data.time(1),data.time(end)]);
    
    legend([l0 l1],{'Reflectivity','Power'},'location','northeast');
    title([{[datestr(data.time(1),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(end),'yyyy-mm-dd HH:MM:SS')]};...
        {'Linear atmosphere power and reflectivity'}]);
    
    s5pos=s5.Position;
    s5.Position=[s5pos(1),s5pos(2),s1pos(3),s5pos(4)];
    
    s2=subplot(5,1,5);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    %plot(data.time,data.ICING_LEVEL./1000,'-k','linewidth',2);
    ylabel('Altitude (km)');
    caxis([-25 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_lwc_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end