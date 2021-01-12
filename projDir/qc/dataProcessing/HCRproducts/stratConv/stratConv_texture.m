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

for aa=2:length(caseStart)
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    %% Get data
    
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    
    data.DBZ = [];
    data.FLAG=[];
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
    
%     %% Find melting layer
%     
%     data.dbzMasked=data.DBZ;
%     data.dbzMasked(data.FLAG>1)=nan;
%     
%     findMelt=f_meltLayer(data,400);
%     oneInds=find(findMelt==1);
%     twoInds=find(findMelt==2);
%     threeInds=find(findMelt==3);
%     
    %% Cloud puzzle
    
    cloudPuzzle=f_cloudPuzzle_radial(data);
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    cloudCount=length(uClouds);
    
    %% Stratiform convective partitioning
    dbzText=nan(size(data.DBZ));
    
    pixRad=350;
    goodDataFrac=0.75;
    
    for jj=1:max(uClouds)
        disp(['Cloud ',num2str(jj),' of ',num2str(max(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=jj)=nan;
        
        % Shrink to good data area
        nonNanCols=find(any(~isnan(dbzIn),1));
        dbzIn=dbzIn(:,nonNanCols);
        
        dbzTextOne=f_reflTexture(dbzIn,pixRad,goodDataFrac);
        dbzTextLarge=nan(size(dbzText));
        dbzTextLarge(:,nonNanCols)=dbzTextOne;
        dbzText(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
    end
    %% Plot strat conv
    
    disp('Plotting ...');
    
    close all
    
    f1 = figure('Position',[200 500 1500 900],'DefaultAxesFontSize',12);
    
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
    surf(data.time,data.asl./1000,dbzText,'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([0 10]);
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    colorbar
    grid on
    title('Reflectivity texture')
    s2pos=s2.Position;
    s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];
    
    s3=subplot(3,1,3);
    colMapIn=jet(cloudCount-1);
    % Make order random
    indsCol=randperm(size(colMapIn,1));
    colMapInds=cat(2,indsCol',colMapIn);
    colMapInds=sortrows(colMapInds);
    colMap=cat(1,[0 0 0],colMapInds(:,2:end));
    
    hold on;
    surf(data.time,data.asl./1000,cloudPuzzle,'edgecolor','none');
    view(2);
    s3.Colormap=colMap;
    ylabel('Altitude (km)');
    ylim([0 ylimUpper]);
    xlim([data.time(1),data.time(end)]);
    grid on
    title('Cloud Puzzle')
    caxis([-0.5 cloudCount-1+0.5])
    colorbar
    s3pos=s3.Position;
    s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stratConv_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end