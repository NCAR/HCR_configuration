% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
dataFreq='10hz';

ylimUpper=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/h/eol/romatsch/hcrCalib/clouds/stratConv/'];
figdir=['/home/romatsch/plots/HCR/stratConv/',project,'/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/stratConv_',project,'.txt'];

%dataDir=HCRdir(project,quality,dataFreq);
dataDir=['/run/media/romatsch/RSF0006/rsf/cloudPuzzle/',project,'/10hz/'];

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
    data.TOPO=[];
    data.CLOUD_PUZZLE=[];
    
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    %frq=ncread(fileList{1},'frequency');
    
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
    %cloudCount=length(uClouds);
    
    %% Calculate reflectivity texture and convectivity
    %stratConvNearSurf=nan(size(data.DBZ));
    
    dbzText=nan(size(data.DBZ));
        
    pixRad=50; % Radius over which texture is calculated in pixels. Default is 50.
    %stratConvThresh=4; % Texture above (below) is convective (stratiform)
    
    for jj=1:length(uClouds)
        disp(['Calculating texture for tethered cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
        dbzIn=data.DBZ;
        dbzIn(data.FLAG>1)=nan;
        dbzIn(cloudPuzzle~=uClouds(jj))=nan;
        
        % Shrink to good data area
        nonNanCols=find(any(~isnan(dbzIn),1));
        dbzIn=dbzIn(:,nonNanCols);
                
        dbzTextOne=f_reflTexture(dbzIn,pixRad);
        
        % Fill into large matrix
        dbzTextLarge=nan(size(dbzText));
        dbzTextLarge(:,nonNanCols)=dbzTextOne;
        dbzText(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
        
        % Prepare asl and topo
%         topoIn=data.TOPO(nonNanCols);
%         aslIn=data.asl(:,nonNanCols);
%         aslIn(isnan(dbzIn))=nan;
%         stratConvNearSurfOne=f_stratConvPart_nearSurf(dbzTextOne,dbzIn,stratConvThresh,dbzThresh,aslIn,topoIn);
%         
%         scLarge=nan(size(data.DBZ));
%         scLarge(:,nonNanCols)=stratConvNearSurfOne;
%         stratConvNearSurf(~isnan(scLarge))=scLarge(~isnan(scLarge));
    end
    
    % Convectivity
    upperLim=8; % Upper limit for convectivity mapping. Texture above that will be set to 1.
    convectivity=1/upperLim.*dbzText;
    
    % Basic classification
    stratMixed=0.3;
    mixedConv=0.5;
    
    classBasic=nan(size(data.DBZ));
    classBasic(convectivity<stratMixed)=1; % Stratiform
    classBasic(convectivity>=stratMixed & convectivity<mixedConv)=2;
    classBasic(convectivity>=mixedConv)=3;
    
%     %% Calculate reflectivity texture for elevated convective echo
%     stratConvTot=nan(size(data.DBZ));
%     
%     dbzText2=nan(size(data.DBZ));
%     
%     pixRad=50; % Radius over which texture is calculated in pixels. Default is 350.
%     dbzThresh=-12; % Reflectivity data below this threshold will not be used in the texture calculation
%     stratConvThresh=3.5; % Texture above (below) is convective (stratiform)
%     
%     for jj=1:length(uClouds)
%         disp(['Calculating texture for elevated cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
%         dbzIn=data.DBZ;
%         dbzIn(data.FLAG>1)=nan;
%         dbzIn(cloudPuzzle~=uClouds(jj))=nan;
%         
%         % Shrink to good data area
%         nonNanCols=find(any(~isnan(dbzIn),1));
%         dbzIn=dbzIn(:,nonNanCols);
%                 
%         dbzTextOne=f_reflTexture(dbzIn,pixRad,dbzThresh);
%         dbzTextLarge=nan(size(dbzText2));
%         dbzTextLarge(:,nonNanCols)=dbzTextOne;
%         dbzText2(~isnan(dbzTextLarge))=dbzTextLarge(~isnan(dbzTextLarge));
%         
%         % Prepare asl and topo
%         topoIn=data.TOPO(nonNanCols);
%         aslIn=data.asl(:,nonNanCols);
%         aslIn(isnan(dbzIn))=nan;
%         stratConvElevOne=f_stratConvPart_elev(dbzTextOne,dbzIn,stratConvThresh,dbzThresh,aslIn,topoIn);
%         
%         scLarge=nan(size(data.DBZ));
%         scLarge(:,nonNanCols)=stratConvElevOne;
%         stratConvTot(~isnan(scLarge))=scLarge(~isnan(scLarge));
%     end
%     
%     stratConvTot(stratConvNearSurf==1)=1;
% 
%     %% Separate isolated and embedded
%     stratConv=nan(size(data.DBZ));
%     
%     for jj=1:length(uClouds)
%         disp(['Stratiform/convective isolated/embedded separation for cloud ',num2str(jj),' of ',num2str(length(uClouds))]);
%                        
%         % StratConv for one cloud
%         stratConvPart=stratConvTot;
%         stratConvPart(cloudPuzzle~=uClouds(jj))=nan;
%         % Shrink to good data area
%         nonNanColsD=find(any(~isnan(stratConvPart),1));
%         stratConvPart=stratConvPart(:,nonNanColsD);
%                 
%         % Divide into sub categories
%         stratConvSub=f_embeddedIsolated(stratConvPart);
%         
%         % Back into large matrix
%         scLarge=nan(size(stratConv));
%         scLarge(:,nonNanColsD)=stratConvSub;
%         stratConv(~isnan(scLarge))=scLarge(~isnan(scLarge));
%     end
%     
%     % Replace strat only areas that are small with strat embedded
%     stratMask=zeros(size(stratConv));
%     stratMask(stratConv==20)=1;
%     stratMask=bwareaopen(stratMask,50000);
%     
%     stratConv(stratConv==20 & stratMask==0)=21;
%     
%     % Replace strat embedded areas that are small with strat only
%     stratMask=zeros(size(stratConv));
%     stratMask(stratConv==21)=1;
%     stratMask=bwareaopen(stratMask,5000);
%     
%     stratConv(stratConv==21 & stratMask==0)=20;
%     
%     stratConv(stratConv==22)=20;
%     
%     %% 1D stratiform convective partitioning
%     stratConv1D=f_stratConv1Dperc(stratConv);
%    
%     %% Small clouds that are 0 in cloud puzzle
%     stratConv(isnan(stratConv) & data.FLAG==1)=30;
%    
%     %% Plot strat conv
%     
%     disp('Plotting ...');
%     
%     close all
%     
%     stratConvPlot=stratConv;
%     stratConvPlot(stratConv==20)=14;
%     stratConvPlot(stratConv==21)=15;
%     stratConvPlot(stratConv==30)=16;
%     
%     colMapSC=[1,0,0;
%         1,0,1;
%         1,0.5,0;
%         0.5,0,1;
%         0,0,1;
%         0,0.7,1;
%         0,0,0];
%
colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    1 0 0];

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
hold on
surf(data.time,data.asl./1000,dbzText,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 12]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
colorbar
grid on
title('Reflectivity texture')
s2pos=s2.Position;
s2.Position=[s2pos(1),s2pos(2),s1pos(3),s2pos(4)];

s3=subplot(4,1,3);

hold on
surf(data.time,data.asl./1000,convectivity,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 1]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
colorbar
grid on
title('Convectivity')
s3pos=s3.Position;
s3.Position=[s3pos(1),s3pos(2),s1pos(3),s3pos(4)];

s4=subplot(4,1,4);

hold on
surf(data.time,data.asl./1000,classBasic,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0.5 3.5]);
colormap(s4,[0,0,1;0.32,0.78,0.59;1 0 0]);
cb=colorbar;
cb.Ticks=[1,2,3];
cb.TickLabels=cat(2,{'Stratiform','Mixed','Convective'});
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
grid on
s4pos=s4.Position;
s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];

linkaxes([s1 s2 s3 s4],'xy');

%
%     s5=subplot(30,1,30);
%
%     hold on
%     scat1=scatter(data.time,ones(size(data.time)),10,stratConv1D,'filled');
%     set(gca,'clim',[0,1]);
%     set(gca,'YTickLabel',[]);
%     %s5.Colormap=;
%     xlim([data.time(1),data.time(end)]);
%     s5pos=s5.Position;
%     s5.Position=[s5pos(1),s5pos(2)-0.023,s1pos(3),s5pos(4)];
%
%     s4=subplot(4,1,4);
%
%     hold on
%     surf(data.time,data.asl./1000,stratConvPlot,'edgecolor','none');
%     view(2);
%     ylabel('Altitude (km)');
%     caxis([0 10]);
%     ylim([0 ylimUpper]);
%     xlim([data.time(1),data.time(end)]);
%     s4.Colormap=colMapSC;
%     caxis([9.5 16.5]);
%     cb=colorbar;
%     cb.Ticks=10:16;
%     cb.TickLabels={'Isolated tethered conv.','Embedded tethered conv.','Isolated elevated conv.',...
%         'Embedded elevated conv.','Stratiform','Strat. with embedded conv.',...
%         'Small'};
%     set(gca,'XTickLabel',[]);
%     grid on
%     title('Stratiform/convective partitioning')
%     s4pos=s4.Position;
%     s4.Position=[s4pos(1),s4pos(2),s1pos(3),s4pos(4)];
%
%     linkaxes([s1 s2 s3 s4],'xy');
%
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_stratConvBasic_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')

end