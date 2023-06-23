% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';
whichModel='era5';

minPixNumUW=5;
largeUW=0; % Set to 1 when we want to use only larges particles

HCRrangePix=5;
HCRtimePix=20;

plotOn=1;
showPlot='off';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

% if strcmp(project,'otrec')
%     indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
% elseif strcmp(project,'socrates')
%     indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
% end

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

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
    'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

cscale_hcr_2=[1 0 0;0 1 0;0 0 1];
units_str_hcr_2={'Liquid','Mixed','Frozen'};

varNames={'numLiqHCR','numIceHCR','numAllHCR','pidHCR','numLiqLargestP','numIceLargestP','numAllLargestP','sizeLargestP'};
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
        
        data.DBZ_MASKED = [];
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
        
        
        %% Find largest
        countAllFlip=flipud(countAll);
        cumSumAll=cumsum(countAllFlip);
        cumSumAllBack=flipud(cumSumAll);
        
        indMat=repmat((1:size(countAll,1))',1,size(countAll,2));
        indMat(countAll==0)=nan;
        indMat(cumSumAllBack<minPixNumUW)=nan;

        if largeUW
            rowsGood=max(indMat,[],1,'omitnan');
        else
            rowsGood=ones(1,size(indMat,2));
            rowsGood(sumAll==0)=nan;
        end

        colsGood=find(~isnan(rowsGood));
        rowsGood=rowsGood(colsGood);
        
        %goodInds=sub2ind(size(countAll),rowsGood,colsGood);
        
        numLiqLargest=nan(1,size(countAll,2));
        numIceLargest=nan(1,size(countAll,2));
        numAllLargest=nan(1,size(countAll,2));
        sizeLargest=nan(1,size(countAll,2));
        
        for jj=1:length(colsGood)
            numLiqLargest(colsGood(jj))=sum(countLiq(rowsGood(jj):end,colsGood(jj)));
            numIceLargest(colsGood(jj))=sum(countIce(rowsGood(jj):end,colsGood(jj)));
            numAllLargest(colsGood(jj))=sum(countAll(rowsGood(jj):end,colsGood(jj)));
            
            addSizes=sizeMM(rowsGood(jj):end).*countAll(rowsGood(jj):end,colsGood(jj));
            sizeLargest(colsGood(jj))=sum(addSizes)./numAllLargest(colsGood(jj));
        end
               
        %% Calculate HCR liquid fraction
        
        hcrLiqIce=nan(size(data.PID));
        hcrLiqIce(data.PID<=6)=1;
        hcrLiqIce(data.PID==7)=2;
        hcrLiqIce(data.PID>=8 & data.PID<10)=3;
        
        liqFrac_HCR_P=nan(length(ptime),6);
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
                    hcrPID=data.PID(18:18+HCRrangePix,hcrIndCols);
                    pidOut=mode(reshape(hcrPID,1,[]));
                    addFrac=[liqNum/allNum,liqFrac(goodIndsP(jj)),liqNum,iceNum,allNum,pidOut];
                    liqFrac_HCR_P((goodIndsP(jj)),:)=addFrac;
                end
            end
        end
        
        outTable=timetable(ptime,liqFrac_HCR_P(:,3),liqFrac_HCR_P(:,4),liqFrac_HCR_P(:,5),liqFrac_HCR_P(:,6),...
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
            surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
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
            s1.Position=[s1pos(1),s1pos(2),s1pos(3),s1pos(4)];
            
            s2=subplot(4,1,2);
            
            hold on
            surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
            view(2);
            colormap(s2,cscale_hcr);
            cb=colorbar;
            cb.Ticks=1:11;
            cb.TickLabels=units_str_hcr;
            ylabel('Altitude (km)');
            title(['HCR particle ID']);
            
            plot(data.time,data.altitude./1000,'-k','linewidth',2);
            
            caxis([.5 11.5]);
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
            plot(ptime,outTable.numLiqLargestP./outTable.numAllLargestP,'-b','linewidth',2);
            plot(ptime,outTable.numLiqHCR./outTable.numAllHCR,'-g','linewidth',2);
            ylim([0 1]);
            yticks(0:0.1:1);
            ylabel('Liquid fraction')
            
            yyaxis right
            set(gca,'YColor','k');
            plot(ptime,sizeLargest,'-k','linewidth',2);
            ylabel('Size (mm)');
            ylim([0 3]);
            yticks(0:0.3:3);
            xlim([data.time(1),data.time(end)]);
            grid on
            title('Liquid fraction and size of largest particles');
            legend('UW','HCR','Size');
            s4pos=s4.Position;
            s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
            
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_pidUW_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
            
        end
     
        startTime=endTime;
    end
end

liqFracHCRall=outTableAll.numLiqHCR./outTableAll.numAllHCR;

liqFracPallL=outTableAll.numLiqLargestP./outTableAll.numAllLargestP;

corrCoeffL=corrcoef(liqFracPallL,liqFracHCRall,'Rows','complete');

%% Hit miss table 1

lowBound=0.1:0.1:0.5;
highBound=fliplr(0.5:0.1:0.9);

xvalues = {'Ice','Mixed','Liquid'};

close all

for ii=1:length(lowBound)
    
    % Hit miss table
    
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

    if ii==3
        save([figdir,'hitMissTable.mat'],'hmTableL');
    end
end

%% Size of different PID classes

pidSize=nan(9,1);
pidStd=nan(9,1);

for ii=1:9
    pidTI=find(outTableAll.pidHCR==ii);
    sizeT=outTableAll.sizeLargestP(pidTI);
    pidSize(ii)=mean(sizeT,'omitnan');
    pidStd(ii)=std(sizeT,'omitnan');
end

close all

f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12,'visible','on','renderer','painters');
hold on

errorbar(1:9,pidSize.*1000,pidStd.*1000,'sk','MarkerSize',1,'linewidth',1.5,'capsize',15);
scatter(1:9,pidSize.*1000,150,cscale_hcr(1:9,:),'filled','MarkerEdgeColor','k','linewidth',2);
%errorbar(1:9,pidSize,pidStd);
xlim([0 10]);
xticks(1:9);
xticklabels(units_str_hcr);
yticks(0:100:4000);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
box on

ylabel('Particle size (\mum)');
title('Mean particle size');

set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_meanPartSizePID.png'],'-dpng','-r0');

save([figdir,'sizes.mat'],'pidSize','pidStd');

%% Save table

save([figdir,'compTable.mat'],'outTableAll');

%% Scatter plot
centers={0.05:0.1:0.95 0.05:0.1:0.95};
cm=jet(230);

close all
f1 = figure('Position',[200 200 1400 1000],'DefaultAxesFontSize',12,'visible','on','renderer','painters');
colormap(cm(30:200,:));

for ii=1:9
    pidTI=find(outTableAll.pidHCR==ii);
    
    lfU=outTableAll.numLiqHCR(pidTI)./outTableAll.numAllHCR(pidTI);
    lfH=outTableAll.numLiqLargestP(pidTI)./outTableAll.numAllLargestP(pidTI);
    
    lfUH=cat(2,lfU,lfH);
    lfUH(any(isnan(lfUH),2),:)=[];
    
    subplot(3,3,ii)
    hist3(lfUH,'Ctrs',centers,'CdataMode','auto','edgecolor','none');
    view(2)
    xlim([0 1])
    ylim([0 1])
    colorbar
    
    title([units_str_hcr{ii},' (',num2str(size(lfUH,1)),')']);
    
    xlabel('HCR')
    ylabel('UWILD')
end

set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_heatMapCats.png'],'-dpng','-r0');

f1 = figure('Position',[200 200 600 500],'DefaultAxesFontSize',12,'visible','on','renderer','painters');
colormap(cm(30:200,:));

liqFracUW_HCR=cat(2,liqFracPallL,liqFracHCRall);
liqFracUW_HCR(any(isnan(liqFracUW_HCR),2),:)=[];

hist3(liqFracUW_HCR,'Ctrs',centers,'CdataMode','auto');
view(2)
xlim([0 1])
ylim([0 1])
colorbar
title(['All (',num2str(size(liqFracUW_HCR,1)),')']);

xlabel('HCR')
ylabel('UWILD')

set(gcf,'PaperPositionMode','auto')
print([figdir,project,'_heatMapAll.png'],'-dpng','-r0');