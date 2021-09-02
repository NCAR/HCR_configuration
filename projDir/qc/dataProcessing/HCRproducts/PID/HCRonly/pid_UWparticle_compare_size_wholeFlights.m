% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
% dataFreq='10hz';
% qcVersion='v2.1';
whichModel='era5';

smallMid=0.14; % Threshold in mm for small to mid particles
midLarge=0.33;

HCRrangePix=10;
HCRtimePix=20;

plotOn=1;
showPlot='off';
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

figdir=[indir(1:end-5),'pidPlots/comparePID_UW_wholeFlights/'];

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid','Mixed Phase','Large Frozen','Small Frozen'};

cscale_hcr_2=[1 0 0;0 1 0;0 0 1];
units_str_hcr_2={'Liquid','Mixed','Frozen'};

varNames={'numLiqHCR','numIceHCR','numAllHCR','numLiqLargestP','numIceLargestP','numAllLargest','sizeLargestP'};
outTableAll=[];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:14
    disp(['Flight ',num2str(aa)]);
    disp('Loading data ...')
    
    %% Get particle data
    
    % Get data
    flightFile=[particleDir,partFiles(aa).name];
    ptimeIn=ncread(flightFile,'time');
    ptime1=partFileTimes(aa)+seconds(ptimeIn);
    
    binEdges=ncread(flightFile,'bin_edges');
    binWidthCM=ncread(flightFile,'bin_width');
    
    sizeMM=(binEdges(1:end-1)+binEdges(2:end))./2;
    
    countLiq1=ncread(flightFile,'count_darea_liq_ml');
    countIce1=ncread(flightFile,'count_darea_ice_ml');
    countAll1=ncread(flightFile,'count_darea_all');
        
    startTime1=datetime(caseList(aa,1:6));
    endTime1=datetime(caseList(aa,7:12));
    
    startTime=startTime1;
    
    while startTime<endTime1
        endTime=startTime+minutes(30);
        
        %% Sub times for particles
        ptimeInds=find(ptime1>=startTime & ptime1<=endTime);
        
        countAll=countAll1(:,ptimeInds);
        
        if sum(sum(countAll))<10
            startTime=endTime;
            continue
        end
        
        ptime=ptime1(ptimeInds);        
        countLiq=countLiq1(:,ptimeInds);
        countIce=countIce1(:,ptimeInds);
                
        sumAll=sum(countAll,1);
          
        liqFrac=sum(countLiq,1)./sumAll;
        
        %% Get HCR data
        
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if isempty(fileList)
            startTime=endTime;
            continue
        end
        
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
        
        if isempty(data.time)
            startTime=endTime;
            continue
        end
        
        data.DBZ(data.FLAG>1)=nan;
        
        %% Find largest
        
        indMat=repmat((1:size(countAll,1))',1,size(countAll,2));
        indMat(countAll==0)=nan;
        
        rowsGood=max(indMat,[],1,'omitnan');
        colsGood=find(~isnan(rowsGood));
        rowsGood=rowsGood(colsGood);
        
        goodInds=sub2ind(size(countAll),rowsGood,colsGood);
        
        numLiqLargest=nan(1,size(countAll,2));
        numIceLargest=nan(1,size(countAll,2));
        numAllLargest=nan(1,size(countAll,2));
        
        numLiqLargest(colsGood)=countLiq(goodInds);
        numIceLargest(colsGood)=countIce(goodInds);
        numAllLargest(colsGood)=countAll(goodInds);
        
        sizeLarge=sizeMM(rowsGood);
        sizeLargest=nan(1,size(countAll,2));
        sizeLargest(colsGood)=sizeLarge;
        
        %% Calculate HCR liquid fraction
        
        hcrLiqIce=nan(size(data.PID));
        hcrLiqIce(data.PID<=6)=1;
        hcrLiqIce(data.PID==7)=2;
        hcrLiqIce(data.PID>=8)=3;
        
        liqFrac_HCR_P=nan(length(ptime),5);
        goodIndsP=find(~isnan(liqFrac));
        
        for jj=1:length(goodIndsP)
            hcrInd=find(data.time==ptime(goodIndsP(jj)));
            if hcrInd>HCRtimePix & hcrInd+HCRtimePix<=length(data.time)
                hcrIndCols=hcrInd-HCRtimePix:hcrInd+HCRtimePix;
                hcrParts=hcrLiqIce(18:18+HCRrangePix,hcrIndCols);
                liqNum=sum(sum(hcrParts==1));
                mixNum=sum(sum(hcrParts==2));
                iceNum=sum(sum(hcrParts==3));
                if mixNum>0
                    liqNum=liqNum+ceil(mixNum/2);
                    iceNum=iceNum+floor(mixNum/2);
                end
                allNum=sum(sum(~isnan(hcrParts)));                
                
                if allNum>50
                    addFrac=[liqNum/allNum,liqFrac(goodIndsP(jj)),liqNum,iceNum,allNum];
                    liqFrac_HCR_P((goodIndsP(jj)),:)=addFrac;
                end
            end
        end
        
        outTable=timetable(ptime,liqFrac_HCR_P(:,3),liqFrac_HCR_P(:,4),liqFrac_HCR_P(:,5),...
            numLiqLargest',numIceLargest',numAllLargest',sizeLargest','VariableNames',varNames);
        
        outTableAll=cat(1,outTableAll,outTable);
        %% Plot 1
        
        if plotOn
            disp('Plotting ...');
            
            ttAlt=timetable(data.time',data.altitude');
            ttP=timetable(ptime);
            ttSync=synchronize(ttP,ttAlt,'first');
            
            liqFracLargest=numLiqLargest./numAllLargest;
            
            cR=[linspace(1,0,51)',zeros(51,2)];
            cB=[zeros(50,2),linspace(0,1,50)'];
            colmapL=flipud(cat(1,cR,cB,[1 1 1]));
            
            liqFracPlotL=liqFracLargest;
            liqFracPlotL=round(liqFracPlotL*100);
            liqFracPlotL=liqFracPlotL+2;
            liqFracPlotL(isnan(liqFracLargest))=1;
            col1DL=colmapL(liqFracPlotL,:);
            
            close all
            
            ylims=[max([0,floor(min(data.altitude./1000))]),ceil(max(data.altitude./1000))];
            
            f1 = figure('Position',[200 500 1800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
            
            s1=subplot(4,1,1);
            
            colormap jet
            
            hold on
            surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([-35 25]);
            ylim(ylims);
            xlim([data.time(1),data.time(end)]);
            colorbar
            grid on
            title('Reflectivity (dBZ)')
            plot(data.time,data.altitude./1000,'-k','linewidth',2);
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
            ylim(ylims);
            xlim([data.time(1),data.time(end)]);
            
            grid on
            box on
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            s3=subplot(4,1,3);
            
            hold on
            surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
            view(2);
            colormap(s3,cscale_hcr_2);
            cb=colorbar;
            cb.Ticks=6:8;
            cb.TickLabels=units_str_hcr_2;
            ylabel('Altitude (km)');
            title(['HCR particle ID']);
            
            scatter(ptime,ttSync.Var1./1000,20,col1DL,'filled');
            set(gca,'clim',[0,1]);
            
            caxis([5.5 8.5]);
            ylim(ylims);
            xlim([data.time(1),data.time(end)]);
            
            grid on
            box on
            s3pos=s3.Position;
            s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
                        
            s4=subplot(4,1,4);
            
            hold on
            plot(ptime,sizeLargest,'-k','linewidth',2);
            ylabel('Size (mm)');
            ylim([0 max(sizeMM)]);
            xlim([data.time(1),data.time(end)]);
            grid on
            title('Size of largest particles');
            s4pos=s4.Position;
            s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_pidUW_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
            
        end
     
        startTime=endTime;
    end
end

liqFracPall=outTableAll.numLiqP./outTableAll.numAllP;
liqFracHCRall=outTableAll.numLiqHCR./outTableAll.numAllHCR;

liqFracPallL=outTableAll.numLiqLargestP./outTableAll.numAllLargest;

corrCoeff=corrcoef(liqFracPall,liqFracHCRall,'Rows','complete');
corrCoeffL=corrcoef(liqFracPallL,liqFracHCRall,'Rows','complete');

%% Hit miss table 1

lowBound=0.1:0.1:0.5;
highBound=fliplr(0.5:0.1:0.9);

close all

for ii=1:length(lowBound)
    
    % HCR first, particles second in name
    iceIce=length(find(liqFracHCRall<lowBound(ii) & liqFracPall<lowBound(ii)));
    iceMix=length(find(liqFracHCRall<lowBound(ii) & liqFracPall>=lowBound(ii) & liqFracPall<=highBound(ii)));
    iceLiq=length(find(liqFracHCRall<lowBound(ii) & liqFracPall>highBound(ii)));
    
    mixIce=length(find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPall<lowBound(ii)));
    mixMix=length(find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPall>=lowBound(ii) & liqFracPall<=highBound(ii)));
    mixLiq=length(find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPall>highBound(ii)));
    
    liqIce=length(find(liqFracHCRall>highBound(ii) & liqFracPall<lowBound(ii)));
    liqMix=length(find(liqFracHCRall>highBound(ii) & liqFracPall>=lowBound(ii) & liqFracPall<=highBound(ii)));
    liqLiq=length(find(liqFracHCRall>highBound(ii) & liqFracPall>lowBound(ii)));
    
    hmTable=[iceIce,mixIce,liqIce;iceMix,mixMix,liqMix;iceLiq,mixLiq,liqLiq];
    
    hmNorm=hmTable./sum(sum(hmTable)).*100;
    
    xvalues = {'Ice','Mixed','Liquid'};
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmNorm);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['All. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Correct: ',...
        num2str(hmNorm(1,1)+hmNorm(2,2)+hmNorm(3,3),3),'%. Correlation: ',num2str(corrCoeff(2,1),2),'.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_All_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
    % Hit miss table 2
    
    % HCR first, particles second in name
    iceIceLinds=find(liqFracHCRall<lowBound(ii) & liqFracPallL<lowBound(ii));
    iceIceL=length(iceIceLinds);
    iceMixLinds=find(liqFracHCRall<lowBound(ii) & liqFracPallL>=lowBound(ii) & liqFracPallL<=highBound(ii));
    iceMixL=length(iceMixLinds);
    iceLiqLinds=find(liqFracHCRall<lowBound(ii) & liqFracPallL>highBound(ii));
    iceLiqL=length(iceLiqLinds);
    
    mixIceLinds=find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPallL<lowBound(ii));
    mixIceL=length(mixIceLinds);
    mixMixLinds=find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPallL>=lowBound(ii) & liqFracPallL<=highBound(ii));
    mixMixL=length(mixMixLinds);
    mixLiqLinds=find(liqFracHCRall>=lowBound(ii) & liqFracHCRall<=highBound(ii) & liqFracPallL>highBound(ii));
    mixLiqL=length(mixLiqLinds);
    
    liqIceLinds=find(liqFracHCRall>highBound(ii) & liqFracPallL<lowBound(ii));
    liqIceL=length(liqIceLinds);
    liqMixLinds=find(liqFracHCRall>highBound(ii) & liqFracPallL>=lowBound(ii) & liqFracPallL<=highBound(ii));
    liqMixL=length(liqMixLinds);
    liqLiqLinds=find(liqFracHCRall>highBound(ii) & liqFracPallL>lowBound(ii));
    liqLiqL=length(liqLiqLinds);
    
    hmTableL=[iceIceL,mixIceL,liqIceL;iceMixL,mixMixL,liqMixL;iceLiqL,mixLiqL,liqLiqL];
    
    hmNormL=hmTableL./sum(sum(hmTableL)).*100;
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmNormL);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW largest particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['Largest. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Correct: ',...
        num2str(hmNormL(1,1)+hmNormL(2,2)+hmNormL(3,3),3),'%. Correlation: ',num2str(corrCoeffL(2,1),2),'.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_Largest_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
    % Temperature
    hmTempL=[mean(outTableAll.tempHCR(iceIceLinds),'omitnan'),mean(outTableAll.tempHCR(mixIceLinds),'omitnan'),mean(outTableAll.tempHCR(liqIceLinds),'omitnan');...
        mean(outTableAll.tempHCR(iceMixLinds),'omitnan'),mean(outTableAll.tempHCR(mixMixLinds),'omitnan'),mean(outTableAll.tempHCR(liqMixLinds),'omitnan');...
        mean(outTableAll.tempHCR(iceLiqLinds),'omitnan'),mean(outTableAll.tempHCR(mixLiqLinds),'omitnan'),mean(outTableAll.tempHCR(liqLiqLinds),'omitnan')];
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmTempL);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW largest particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['Largest. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Mean HCR temperature.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_Largest_Temp_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
    % Reflectivity
    hmReflL=[mean(outTableAll.reflHCR(iceIceLinds),'omitnan'),mean(outTableAll.reflHCR(mixIceLinds),'omitnan'),mean(outTableAll.reflHCR(liqIceLinds),'omitnan');...
        mean(outTableAll.reflHCR(iceMixLinds),'omitnan'),mean(outTableAll.reflHCR(mixMixLinds),'omitnan'),mean(outTableAll.reflHCR(liqMixLinds),'omitnan');...
        mean(outTableAll.reflHCR(iceLiqLinds),'omitnan'),mean(outTableAll.reflHCR(mixLiqLinds),'omitnan'),mean(outTableAll.reflHCR(liqLiqLinds),'omitnan')];
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmReflL);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW largest particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['Largest. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Mean HCR reflectivity.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_Largest_Refl_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
    % Std Reflectivity
    hmStdReflL=[mean(outTableAll.stdReflHCR(iceIceLinds),'omitnan'),mean(outTableAll.stdReflHCR(mixIceLinds),'omitnan'),mean(outTableAll.stdReflHCR(liqIceLinds),'omitnan');...
        mean(outTableAll.stdReflHCR(iceMixLinds),'omitnan'),mean(outTableAll.stdReflHCR(mixMixLinds),'omitnan'),mean(outTableAll.stdReflHCR(liqMixLinds),'omitnan');...
        mean(outTableAll.stdReflHCR(iceLiqLinds),'omitnan'),mean(outTableAll.stdReflHCR(mixLiqLinds),'omitnan'),mean(outTableAll.stdReflHCR(liqLiqLinds),'omitnan')];
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmStdReflL);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW largest particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['Largest. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Mean HCR stddev of reflectivity.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_Largest_ReflStd_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
    % Gradient Reflectivity
    hmGradReflL=[mean(outTableAll.gradReflHCR(iceIceLinds),'omitnan'),mean(outTableAll.gradReflHCR(mixIceLinds),'omitnan'),mean(outTableAll.gradReflHCR(liqIceLinds),'omitnan');...
        mean(outTableAll.gradReflHCR(iceMixLinds),'omitnan'),mean(outTableAll.gradReflHCR(mixMixLinds),'omitnan'),mean(outTableAll.gradReflHCR(liqMixLinds),'omitnan');...
        mean(outTableAll.gradReflHCR(iceLiqLinds),'omitnan'),mean(outTableAll.gradReflHCR(mixLiqLinds),'omitnan'),mean(outTableAll.gradReflHCR(liqLiqLinds),'omitnan')];
    
    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on');
    
    h=heatmap(xvalues,xvalues,hmGradReflL);
    ax = gca;
    axp = struct(ax);       %you will get a warning
    axp.Axes.XAxisLocation = 'top';
    h.ColorbarVisible = 'off';
    
    h.XLabel = 'HCR';
    h.YLabel = 'UW largest particles';
    h.CellLabelFormat = '%.1f';
    h.Title = ['Largest. Boundaries: ',num2str(lowBound(ii)),', ',num2str(highBound(ii)),'. Mean HCR gradient of reflectivity.'];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stats_Largest_ReflGrad_point',num2str(lowBound(ii)*10),'point',num2str(highBound(ii)*10),'.png'],'-dpng','-r0')
    
end