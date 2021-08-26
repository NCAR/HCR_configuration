% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
% dataFreq='10hz';
% qcVersion='v2.1';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=4;
end

smallMid=0.14; % Threshold in mm for small to mid particles
midLarge=0.33;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

if strcmp(project,'otrec')
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
end

%% Get times of UW data

particleDir='/scr/snow2/rsfdata/projects/socrates/microphysics/UW_IceLiquid/';

partFiles=dir([particleDir,'UW_particle_classifications.1hz.*.nc']);

partFileTimes=[];
for ii=1:length(partFiles)
    thisFile=[particleDir,partFiles(ii).name];
    
    fileInfo=ncinfo(thisFile);
    fileTimeStr=fileInfo.Variables(1).Attributes.Value;
    fileTime=datetime(str2num(fileTimeStr(15:18)),str2num(fileTimeStr(20:21)),str2num(fileTimeStr(23:24)),...
        str2num(fileTimeStr(26:27)),str2num(fileTimeStr(29:30)),str2num(fileTimeStr(32:33)));
    
    partFileTimes=cat(1,partFileTimes,fileTime);
end

%% HCR data

figdir=[indir(1:end-5),'pidPlots/comparePID_UW/'];

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid','Mixed Phase','Large Frozen','Small Frozen'};

cscale_hcr_2=[1 0 0;0 1 0;0 0 1];
units_str_hcr_2={'Liquid','Mixed','Frozen'};

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/pid_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Get HCR data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.FLAG=[];
    data.PID=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
    
    dataVars=dataVars(~cellfun('isempty',dataVars));
    
    %% Get particle data
    
    % Find right flight
    flightInd=max(find(partFileTimes<startTime));
    
    % Get data
    flightFile=[particleDir,partFiles(flightInd).name];
    ptimeIn=ncread(flightFile,'time');
    ptime=partFileTimes(flightInd)+seconds(ptimeIn);
    
    ptimeInds=find(ptime>=startTime & ptime<=endTime);
    ptime=ptime(ptimeInds);
    
    binEdges=ncread(flightFile,'bin_edges');
    binWidthCM=ncread(flightFile,'bin_width');
    
    countLiq=ncread(flightFile,'count_darea_liq_ml');
    countLiq=countLiq(:,ptimeInds);
    countIce=ncread(flightFile,'count_darea_ice_ml');
    countIce=countIce(:,ptimeInds);
    countAll=ncread(flightFile,'count_darea_all');
    countAll=countAll(:,ptimeInds);
    phaseFlip=ncread(flightFile,'phase_flip_counts');
    phaseFlip=phaseFlip(ptimeInds);
    
    %% Divide in small, medium, and large particles
    
    smallInd=find(binEdges<smallMid);
    midInd=find(binEdges>=smallMid & binEdges<midLarge);
    largeInd=find(binEdges>=midLarge);
    largeInd=largeInd(1:end-1);
    
    smallLiq=sum(countLiq(smallInd,:),1);
    midLiq=sum(countLiq(midInd,:),1);
    largeLiq=sum(countLiq(largeInd,:),1);
    
    smallIce=sum(countIce(smallInd,:),1);
    midIce=sum(countIce(midInd,:),1);
    largeIce=sum(countIce(largeInd,:),1);
    
    smallAll=sum(countAll(smallInd,:),1);
    midAll=sum(countAll(midInd,:),1);
    largeAll=sum(countAll(largeInd,:),1);
    
    % Fractions
    smallLiqFrac=smallLiq./smallAll;
    midLiqFrac=midLiq./midAll;
    largeLiqFrac=largeLiq./largeAll;
    
    smallIceFrac=smallIce./smallAll;
    midIceFrac=midIce./midAll;
    largeIceFrac=largeIce./largeAll;
    
    sumAll=sum(countAll,1);
    sumAll(sumAll<5)=nan;
    
    liqFrac=sum(countLiq,1)./sumAll;
    
    %% Plot 1
    
    disp('Plotting ...');
    
    close all
    
    f1 = figure('Position',[200 500 1800 1200],'DefaultAxesFontSize',12,'visible','off');
    
    s1=subplot(4,1,1);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-35 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s1pos=s1.Position;
    
    s2=subplot(4,1,2);
    
    hold on
    surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
    view(2);
    colormap(s2,cscale_hcr);
    cb=colorbar;
    cb.Ticks=1:9;
    cb.TickLabels=units_str_hcr;
    ylabel('Altitude (km)');
    title(['HCR particle ID']);
    
    plot(data.time,data.altitude./1000,'-k','linewidth',2);
    
    caxis([.5 9.5]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    
    grid on
    box on
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(4,1,3);
    
    hold on
    plot(ptime,smallLiq,'-g','linewidth',2);
    plot(ptime,midLiq,'-c','linewidth',2);
    plot(ptime,largeLiq,'-b','linewidth',2);
    plot(ptime,smallIce,'-y','linewidth',2);
    plot(ptime,midIce,'-m','linewidth',2);
    plot(ptime,largeIce,'-r','linewidth',2);
    ylabel('Count');
    %ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    grid on
    title('Particle counts');
    legend('Small liquid','Medium liquid','Large liquid','Small ice','Medium ice','Large ice');
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

    s4=subplot(4,1,4);
    
    hold on
    plot(ptime,smallLiq./sumAll,'-g','linewidth',2);
    plot(ptime,midLiq./sumAll,'-c','linewidth',2);
    plot(ptime,largeLiq./sumAll,'-b','linewidth',2);
    plot(ptime,smallIce./sumAll,'-y','linewidth',2);
    plot(ptime,midIce./sumAll,'-m','linewidth',2);
    plot(ptime,largeIce./sumAll,'-r','linewidth',2);
    ylabel('Fraction');
    ylim([0 1]);
    xlim([data.time(1),data.time(end)]);
    grid on
    title('Particle fractions');
    legend('Small liquid','Medium liquid','Large liquid','Small ice','Medium ice','Large ice');
    s4pos=s4.Position;
    s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_pidUW_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    
    %% Plot 2
    
    cR=[linspace(1,0,50)',zeros(50,2)];
    cB=[zeros(50,2),linspace(0,1,50)'];
    colmapL=flipud(cat(1,cR,cB,[1 1 1]));
    
    liqFracPlot=liqFrac;
    liqFracPlot(isnan(liqFrac))=0;
    liqFracPlot=round(liqFracPlot*100);
    col1D=colmapL(liqFracPlot+1,:);
    
    ttAlt=timetable(data.time',data.altitude');
    ttP=timetable(ptime);
    ttSync=synchronize(ttP,ttAlt,'first');
        
    disp('Plotting 2 ...');
   
    f1 = figure('Position',[200 500 1800 1200],'DefaultAxesFontSize',12,'visible','off');
    
    s1=subplot(3,1,1);
    
    colormap jet
    
    hold on
    surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-35 25]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')
    s1pos=s1.Position;
    
    s2=subplot(3,1,2);
    
    hold on
    surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
    view(2);
    colormap(s2,cscale_hcr_2);
    cb=colorbar;
    cb.Ticks=6:8;
    cb.TickLabels=units_str_hcr_2;
    ylabel('Altitude (km)');
    title(['HCR particle ID']);
    
    scatter(ptime,ttSync.Var1./1000,20,col1D,'filled');
    set(gca,'clim',[0,1]);
    
    caxis([5.5 8.5]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    
    grid on
    box on
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(3,1,3);
    
    hold on
    plot(ptime,liqFrac,'-r','linewidth',2);
    ylabel('Fraction');
    ylim([0 1]);
    xlim([data.time(1),data.time(end)]);
    grid on
    title('Liquid fraction');
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_pidUW_liquidIce_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end