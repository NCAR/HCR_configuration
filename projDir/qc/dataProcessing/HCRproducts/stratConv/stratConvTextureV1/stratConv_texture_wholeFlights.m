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
    data.MELTING_LAYER=[];
    data.TOPO=[];
    data.ICING_LEVEL=[];
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
    cloudCount=length(uClouds);
    
    %% Calculate reflectivity texture
    dbzText=nan(size(data.DBZ));
    
    pixRad=350; % Radius over which texture is calculated in pixels
    dbzThresh=-10; % Reflectivity data below this threshold will not be used in the texture calculation
    
    for jj=1:length(uClouds)
        disp(['Calculating texture for cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
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
    end
    
    %% Stratiform convective partitioning
    
    % Big clouds
    stratConvThresh=4; % Texture above (below) is convective (stratiform)
    data.MELTING_LAYER(data.MELTING_LAYER<20)=10;
    data.MELTING_LAYER(data.MELTING_LAYER>19)=20;
    
    stratConv=nan(size(dbzText));
    for jj=1:length(uClouds)
        disp(['Stratiform/convective partitioning for cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        
        % Reflectivity for one cloud
        dbzPart=data.DBZ;
        dbzPart(data.FLAG>1)=nan;
        dbzPart(cloudPuzzle~=uClouds(jj))=nan;
        % Shrink to good data area
        nonNanColsD=find(any(~isnan(dbzPart),1));
        dbzPart=dbzPart(:,nonNanColsD);
        
        % Texture for one cloud
        textPart=dbzText;
        textPart(cloudPuzzle~=uClouds(jj))=nan;
        % Shrink to good data area
        textPart=textPart(:,nonNanColsD);
        
        stratConvFun=f_stratConvTexture(textPart,dbzPart,stratConvThresh,dbzThresh);
        
        % Divide into sub categories
        stratConvSub=f_stratConvSub(stratConvFun,data.MELTING_LAYER(:,nonNanColsD));
        
        % Back into large matrix
        scLarge=nan(size(stratConv));
        scLarge(:,nonNanColsD)=stratConvSub;
        stratConv(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    %% 1D stratiform convective partitioning
    stratConv1D=f_stratConv1D(stratConv,data.MELTING_LAYER,data.ICING_LEVEL,data.asl,data.TOPO);
    
    %% Small clouds that are 0 in cloud puzzle
    dbzSmallClouds=data.DBZ;
    dbzSmallClouds(data.FLAG>1)=nan;
    
    [stratConv,stratConv1D]=f_stratConvSmallClouds(stratConv,stratConv1D,cloudPuzzle,dbzSmallClouds,data.MELTING_LAYER);
    
    %% Plot in hourly increments
    
    disp('Plotting ...');
    
    colMapSC=[1,0,0;
        1,0,1;
        0.5,0,1;
        1,0.5,0;
        1,1,0;
        0,0,0.5;
        0,0,1;
        0,0.7,1;
        0,1,1];
    
    startPlot=startTime;
    
    while startPlot<endTime
        endPlot=startPlot+minutes(15);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);
        
        timePlot=data.time(timeInds);
        dbzPlot=data.DBZ(:,timeInds);
        aslPlot=data.asl(:,timeInds);
        puzzlePlot=cloudPuzzle(:,timeInds);
        textPlot=dbzText(:,timeInds);
        
        close all
        
        stratConvPlot=stratConv(:,timeInds);
        stratConvPlot(stratConv(:,timeInds)==20)=15;
        stratConvPlot(stratConv(:,timeInds)==21)=16;
        stratConvPlot(stratConv(:,timeInds)==22)=17;
        stratConvPlot(stratConv(:,timeInds)==23)=18;
        
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
            minCloud=uClouds(2);
            maxCloud=uClouds(end);
            
            puzzlePlot(puzzlePlot==0)=minCloud-1;
                        
            colMapIn=jet(maxCloud-minCloud+1);
            % Make order random
            indsCol=randperm(size(colMapIn,1));
            colMapInds=cat(2,indsCol',colMapIn);
            colMapInds=sortrows(colMapInds);
            colMap=cat(1,[0 0 0],colMapInds(:,2:end));
            
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
            
            hold on;
            surf(timePlot,aslPlot./1000,puzzlePlot,'edgecolor','none');
            view(2);
            s2.Colormap=colMap;
            ylabel('Altitude (km)');
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            grid on
            title('Cloud Puzzle')
            caxis([minCloud-1.5 maxCloud+0.5])
            cb=colorbar;
            cb.TickLabels={''};
            s2pos=s2.Position;
            s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
            
            s3=subplot(4,1,3);
            
            hold on
            surf(timePlot,aslPlot./1000,textPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([0 10]);
            ylim([0 ylimUpper]);
            xlim([timePlot(1),timePlot(end)]);
            colorbar
            grid on
            title('Reflectivity texture')
            s3pos=s3.Position;
            s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
            
            s5=subplot(30,1,30);
            timeConv=timePlot(stratConv1D(timeInds)==1);
            conv1D=ones(1,length(timeConv));
            timeStrat=timePlot(stratConv1D(timeInds)==2);
            strat1D=ones(1,length(timeStrat));
            
            hold on
            scatter(timeStrat,strat1D,10,'b','filled');
            scatter(timeConv,conv1D,10,'r','filled');
            set(gca,'YTickLabel',[]);
            
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
            caxis([9.5 18.5]);
            cb=colorbar;
            cb.Ticks=10:18;
            cb.TickLabels={'Isolated conv.','Warm embedded conv.','Cold embedded conv.',...
                'Warm small conv.','Cold small conv.',...
                'Isolated strat.','Strat. with embedded conv.',...
                'Cold small strat.','Warm small strat.'};
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