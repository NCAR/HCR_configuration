% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

ylimUpper=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/h/eol/romatsch/hcrCalib/clouds/stratConv/'];
figdir=['/home/romatsch/plots/HCR/stratConv/texture/',project,'/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/stratConv_',project,'.txt'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

% Loop through cases

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
    data.FLAG=[];
    data.ICING_LEVEL=[];
    %data.MELTING_LAYER=[];
    
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
    disp('Making cloud puzzle');
    
    cloudPuzzle=f_cloudPuzzle_radial(data);
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    cloudCount=length(uClouds);
    
    %% Calculate reflectivity texture
    dbzText=nan(size(data.DBZ));
    
    pixRad=350; % Radius over which texture is calculated in pixels
    dbzThresh=-10; % Reflectivity data below this threshold will not be used in the texture calculation
    
    for jj=1:max(uClouds)
        disp(['Calculating texture for cloud ',num2str(jj),' of ',num2str(max(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=jj)=nan;
        
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
    
    stratConv=nan(size(dbzText));
    for jj=1:max(uClouds)
        disp(['Stratiform/convective partitioning for cloud ',num2str(jj),' of ',num2str(max(uClouds))]);
                    
        % Reflectivity for one cloud
        dbzPart=data.DBZ;
        dbzPart(data.FLAG>1)=nan;
        dbzPart(cloudPuzzle~=jj)=nan;
        % Shrink to good data area
        nonNanColsD=find(any(~isnan(dbzPart),1));
        dbzPart=dbzPart(:,nonNanColsD);
        
        % Texture for one cloud
        textPart=dbzText;
        textPart(cloudPuzzle~=jj)=nan;
        % Shrink to good data area
        textPart=textPart(:,nonNanColsD);
        
        stratConvFun=f_stratConvTexture(textPart,dbzPart,stratConvThresh,dbzThresh);
        
        % Back into large matrix
        scLarge=nan(size(stratConv));
        scLarge(:,nonNanColsD)=stratConvFun;
        stratConv(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    %% 1D stratiform convective partitioning
   
    stratConv1D=f_stratConv1D(stratConv,data.asl,data.ICING_LEVEL);
   
    %% Small clouds that are 0 in cloud puzzle
    dbzSmallClouds=data.DBZ;
    dbzSmallClouds(data.FLAG>1)=nan;
    
    [stratConv,stratConv1D]=f_stratConvSmallClouds(stratConv,stratConv1D,cloudPuzzle,dbzSmallClouds,data.asl,data.ICING_LEVEL);
    %% Plot strat conv
    
    disp('Plotting ...');
    
    close all
    
    f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
    
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
    colMapIn=jet(cloudCount-1);
    % Make order random
    indsCol=randperm(size(colMapIn,1));
    colMapInds=cat(2,indsCol',colMapIn);
    colMapInds=sortrows(colMapInds);
    colMap=cat(1,[0 0 0],colMapInds(:,2:end));
    
    hold on;
    surf(data.time,data.asl./1000,cloudPuzzle,'edgecolor','none');
    view(2);
    s2.Colormap=colMap;
    ylabel('Altitude (km)');
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    grid on
    title('Cloud Puzzle')
    caxis([-0.5 cloudCount-1+0.5])
    colorbar
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(4,1,3);
        
    hold on
    surf(data.time,data.asl./1000,dbzText,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity texture')
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    s5=subplot(30,1,30);
    timeConv=data.time(stratConv1D==1);
    conv1D=ones(1,length(timeConv));
    timeStrat=data.time(stratConv1D==2);
    strat1D=ones(1,length(timeStrat));
    
    hold on
    scatter(timeStrat,strat1D,10,'b','filled');
    scatter(timeConv,conv1D,10,'r','filled');
    set(gca,'YTickLabel',[]);
    
    xlim([data.time(1),data.time(end)]);
    s5pos=s5.Position;
    s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];
    
    s4=subplot(4,1,4);
        
    hold on
    surf(data.time,data.asl./1000,stratConv,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    s4.Colormap=[1 0 0;0 0 1];
    caxis([0.5 2.5]);
    cb=colorbar;
    cb.Ticks=[1 2];
    cb.TickLabels={'Convective','Stratiform'};
    set(gca,'XTickLabel',[]);
    grid on
    title('Stratiform/convective partitioning')
    s4pos=s4.Position;
    s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
        
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stratConv_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end