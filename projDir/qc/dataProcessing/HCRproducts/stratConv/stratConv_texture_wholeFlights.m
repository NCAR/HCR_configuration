% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';
whichModel='era5';

if strcmp(project,'otrec')
    ylimUpper=15;
else
    ylimUpper=10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/cloudPuzzle/',project,'/10hz/'];

outdir=['/run/media/romatsch/RSF0006/rsf/stratConvHCR/',project,'Mat/'];

%figdir=['/h/eol/romatsch/hcrCalib/clouds/stratConv/'];
figdir=['/home/romatsch/plots/HCR/stratConv/',project,'/plotAll/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=4:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.FLAG=[];
    data.TOPO=[];
    data.CLOUD_PUZZLE=[];
    
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
    
    %% Cloud puzzle
    %disp('Making cloud puzzle');
    
    %cloudPuzzle=f_cloudPuzzle_radial(data);
    cloudPuzzle=data.CLOUD_PUZZLE;
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    uClouds(uClouds==0)=[];
        
    %% Calculate reflectivity texture for near surface convective echo
    stratConvNearSurf=nan(size(data.DBZ));
    
    dbzText=nan(size(data.DBZ));
    
    pixRad=200; % Radius over which texture is calculated in pixels. Default is 350.
    dbzThresh=-10; % Reflectivity data below this threshold will not be used in the texture calculation
    stratConvThresh=4; % Texture above (below) is convective (stratiform)
    
    for jj=1:length(uClouds)
        disp(['Calculating texture for tethered cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=uClouds(jj))=nan;
        
        % Shrink to good data area
        nonNanCols=find(any(~isnan(dbzIn),1));
        dbzIn=dbzIn(:,nonNanCols);
        
        dbzTextOne=f_reflTexture(dbzIn,pixRad,dbzThresh);
        dbzTextLarge=nan(size(dbzText));
        dbzTextLarge(:,nonNanCols)=dbzTextOne;
        dbzText(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
        
        % Prepare asl and topo
        topoIn=data.TOPO(nonNanCols);
        aslIn=data.asl(:,nonNanCols);
        aslIn(isnan(dbzIn))=nan;
        stratConvNearSurfOne=f_stratConvPart_nearSurf(dbzTextOne,dbzIn,stratConvThresh,dbzThresh,aslIn,topoIn);
        
        scLarge=nan(size(data.DBZ));
        scLarge(:,nonNanCols)=stratConvNearSurfOne;
        stratConvNearSurf(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    %% Calculate reflectivity texture for elevated convective echo
    stratConvTot=nan(size(data.DBZ));
    
    dbzText2=nan(size(data.DBZ));
    
    pixRad=50; % Radius over which texture is calculated in pixels. Default is 350.
    dbzThresh=-12; % Reflectivity data below this threshold will not be used in the texture calculation
    stratConvThresh=3.5; % Texture above (below) is convective (stratiform)
    
    for jj=1:length(uClouds)
        disp(['Calculating texture for elevated cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=uClouds(jj))=nan;
        
        % Shrink to good data area
        nonNanCols=find(any(~isnan(dbzIn),1));
        dbzIn=dbzIn(:,nonNanCols);
        
        dbzTextOne=f_reflTexture(dbzIn,pixRad,dbzThresh);
        dbzTextLarge=nan(size(dbzText2));
        dbzTextLarge(:,nonNanCols)=dbzTextOne;
        dbzText2(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
        
        % Prepare asl and topo
        topoIn=data.TOPO(nonNanCols);
        aslIn=data.asl(:,nonNanCols);
        aslIn(isnan(dbzIn))=nan;
        stratConvElevOne=f_stratConvPart_elev(dbzTextOne,dbzIn,stratConvThresh,dbzThresh,aslIn,topoIn);
        
        scLarge=nan(size(data.DBZ));
        scLarge(:,nonNanCols)=stratConvElevOne;
        stratConvTot(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    stratConvTot(stratConvNearSurf==1)=1;
    
    %% Separate isolated and embedded
    stratConv=nan(size(data.DBZ));
    
    for jj=1:length(uClouds)
        disp(['Stratiform/convective isolated/embedded separation for cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        
        % StratConv for one cloud
        stratConvPart=stratConvTot;
        stratConvPart(cloudPuzzle~=uClouds(jj))=nan;
        % Shrink to good data area
        nonNanColsD=find(any(~isnan(stratConvPart),1));
        stratConvPart=stratConvPart(:,nonNanColsD);
        
        % Divide into sub categories
        stratConvSub=f_embeddedIsolated(stratConvPart);
        
        % Back into large matrix
        scLarge=nan(size(stratConv));
        scLarge(:,nonNanColsD)=stratConvSub;
        stratConv(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    %% 1D stratiform convective partitioning
    stratConv1D=f_stratConv1Dperc(stratConv);
    
    %% Small clouds that are 0 in cloud puzzle
    stratConv(isnan(stratConv) & data.FLAG==1)=30;
    %% Plot in hourly increments
    
    disp('Plotting ...');
   
    startPlot=startTime;
    
    colMapSC=[1,0,0;
        1,0,1;
        1,0.5,0;
        0.5,0,1;
        0,0,1;
        0,0.7,1;
        0,0,0];
    
    while startPlot<endTime
        
        close all
        
        endPlot=startPlot+minutes(15);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        aslPlot=data.asl(:,timeInds);
        puzzlePlot=cloudPuzzle(:,timeInds);
        textPlot=dbzText(:,timeInds);
        textPlot2=dbzText2(:,timeInds);
        sc1D=stratConv1D(timeInds);
        
        uClouds=unique(puzzlePlot(~isnan(puzzlePlot)));
        if isempty(uClouds)
            startPlot=endPlot;
            continue
        end
        
        if uClouds(1)~=0
            uClouds=[0;uClouds];
        end
        cloudCount=length(uClouds);
        
        if cloudCount>0 & length(uClouds)~=1
            stratConvPlot=stratConv(:,timeInds);
            stratConvPlot(stratConvPlot==20)=14;
            stratConvPlot(stratConvPlot==21)=15;
            stratConvPlot(stratConvPlot==30)=16;
            
            
            f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
            
            s1=subplot(4,1,1);
            
            colormap jet
            
            hold on
            surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([-35 25]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            colorbar
            grid on
            title('Reflectivity (dBZ)')
            s1pos=s1.Position;
            
            s2=subplot(4,1,2);
            hold on
            surf(timePlot,aslPlot./1000,textPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([0 8]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            colorbar
            grid on
            title('Reflectivity texture')
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            s3=subplot(4,1,3);
            
            hold on
            surf(timePlot,aslPlot./1000,textPlot2,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([-2 6]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            colorbar
            grid on
            title('Reflectivity texture')
            s3pos=s3.Position;
            s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
            
            s5=subplot(30,1,30);
            
            hold on
            scat1=scatter(timePlot,ones(size(timePlot)),10,sc1D,'filled');
            set(gca,'clim',[0,1]);
            set(gca,'YTickLabel',[]);
            %s5.Colormap=;
            xlim([timePlot(1),timePlot(end)]);
            s5pos=s5.Position;
            s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];
            
            s4=subplot(4,1,4);
            
            hold on
            surf(timePlot,aslPlot./1000,stratConvPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([0 10]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            s4.Colormap=colMapSC;
            caxis([9.5 16.5]);
            cb=colorbar;
            cb.Ticks=10:16;
            cb.TickLabels={'Isolated tethered conv.','Embedded tethered conv.','Isolated elevated conv.',...
                'Embedded elevated conv.','Stratiform','Strat. with embedded conv.',...
                'Small'};
            set(gca,'XTickLabel',[]);
            grid on
            title('Stratiform/convective partitioning')
            s4pos=s4.Position;
            s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
                        
            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_Flight',num2str(aa),'_stratConv_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        end
        startPlot=endPlot;
    end
    
    %% Save
    disp('Saving stratConv field ...')
    
    save([outdir,whichModel,'.stratConv.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'stratConv');
    
    save([outdir,whichModel,'.stratConv1D.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'stratConv1D');
    
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
end