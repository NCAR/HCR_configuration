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
%     %% Cloud puzzle
%     
%     cloudPuzzle=f_cloudPuzzle_radial(data);
    
    %% Stratiform convective partitioning
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
    
    pixRad=120;
    goodDataFrac=0.75;
    dbzText=f_reflTexture(data,pixRad,goodDataFrac);
    
    %% Plot strat conv
    
    disp('Plotting ...');
    
    %close all
    
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
    
    colmap3=lines;
    colmap3=cat(1,[0 0 0],colmap3);
    
%     hold on
%     surf(data.time,data.asl./1000,stratConv,'edgecolor','none');
%     view(2);
%     colormap(s3,colmap3)
%     ylabel('Altitude (km)');
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     grid on
%     title('Stratiform/convective')
%     s3pos=s3.Position;
%     s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stratConv_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
end