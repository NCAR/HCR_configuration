% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.2';
freqData='10hz'; % 10hz, 100hz, or 2hz

% Determines plot zoom.
if strcmp(project,'otrec')
    ylimits=[-0.3 14];
elseif strcmp(project,'socrates')
    ylimits=[-0.2 6];
elseif strcmp(project,'cset')
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sci/romatsch/meltLayerHCR/',project,'/cases/'];
%figdir=['/home/romatsch/plots/HCR/meltingLayer/selected/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/meltLayer_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);
%indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];


startTime=datetime(2019,9,25,13,41,0);
endTime=datetime(2019,9,25,14,6,0);

%% Load data

disp('Loading data ...');

data=[];

data.DBZ=[];
data.RH=[];
data.FLAG=[];
data.MELTING_LAYER=[];
data.TOPO=[];

dataVars=fieldnames(data);

% Make list of files within the specified time frame
fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if length(fileList)==0
    disp('No data files found.');
    return
end

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for ii=1:length(dataVars)
    if ~isfield(data,dataVars{ii})
        dataVars{ii}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

%% Find melting layer

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);

 %% Cloud puzzle
    %disp('Making cloud puzzle');
    
    cloudPuzzle=f_cloudPuzzle_radial(data);
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    uClouds(uClouds==0)=[];
    cloudCount=length(uClouds);
    
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
    
    % Calculate reflectivity texture for elevated convective echo
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

    % Separate isolated and embedded
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
    
    % Replace strat only areas that are small with strat embedded
    stratMask=zeros(size(stratConv));
    stratMask(stratConv==20)=1;
    stratMask=bwareaopen(stratMask,50000);
    
    stratConv(stratConv==20 & stratMask==0)=21;
    
    % Replace strat embedded areas that are small with strat only
    stratMask=zeros(size(stratConv));
    stratMask(stratConv==21)=1;
    stratMask=bwareaopen(stratMask,5000);
    
    stratConv(stratConv==21 & stratMask==0)=20;
    
    stratConv(stratConv==22)=20;
    
    % Small clouds that are 0 in cloud puzzle
    stratConv(isnan(stratConv) & data.FLAG==1)=30;

%% Plot

timeMat=repmat(data.time,size(data.DBZ,1),1);
dbzMasked=data.DBZ;
dbzMasked(data.FLAG>1)=nan;
humMasked=data.RH;
flagMasked=data.FLAG;

close all

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end

% Resample for plotting
newDBZ=dbzMasked(:,newInds);
newHum=humMasked(:,newInds);
newFLAG=flagMasked(:,newInds);
newASL=data.asl(:,newInds);
newTime=data.time(newInds);

stratConvPlot=stratConv;
stratConvPlot(stratConv==20)=14;
stratConvPlot(stratConv==21)=15;
stratConvPlot(stratConv==30)=16;

newStratConv=stratConvPlot(:,newInds);

colMapSC=[1,0,0;
    1,0,1;
    1,0.5,0;
    0.5,0,1;
    0,0,1;
    0,0.7,1;
    0,0,0];

ytickLabels={'Cloud (1)';'Speckle (2)';'Extinct (3)';'Backlobe (4)';'Out of range (5)';...
    'Bang (6)';'Water (7)';'Land (8)';'Below surf. (9)';...
    'NS cal (10)';'Ant. trans. (11)';'Missing (12)'};

colMask=[0.4,0.8,1;
    0,0,0;
    0.5,0,0.5;
    0,1,0;
    0.2,0.6,0.2;
    1,0,0;
    0,0,0.6;
    0.7065,0.4415,0.2812;
    0.5,0.5,0.5;    
    0.9290,0.8940,0.1250;
    1,0,1;    
    1,0.6,0];

fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1200,920]);

ax1=subplot(4,1,1);
hold on;
sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
view(2);
sub1=colMapDBZ(sub1);
scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');

scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');
scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');
ax = gca;
ax.SortMethod = 'childorder';
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity and melting layer')
grid on
set(gca,'xticklabel',[])
ax1.Position=[0.06 0.765 0.8 0.21];

ax2=subplot(4,1,2);
hold on;
sub1=surf(newTime,newASL./1000,newHum,'edgecolor','none');
caxis([0 100])
ax2.Colormap=flipud(jet);
cb2=colorbar;
set(get(cb2,'Title'),'String','%');
view(2);

ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('ERA5 relative humidity')
grid on
set(gca,'xticklabel',[])
ax2.Position=[0.06 0.525 0.8 0.21];

%%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3=subplot(4,1,3);
hold on;
sub3=surf(newTime,newASL./1000,newFLAG,'edgecolor','none');
view(2);
caxis([1 12]);
colormap(ax3,colMask);
hcb=colorbar;
set(hcb,'ytick',[1.5:0.91:12.5]);
set(hcb,'YTickLabel',ytickLabels);
ylim(ylimits);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Quality flag')
grid on
set(gca,'xticklabel',[])
ax3.Position=[0.06 0.287 0.8 0.21];

ax4=subplot(4,1,4);
s4=surf(newTime,newASL./1000,newStratConv,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([0 10]);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
ax4.Colormap=colMapSC;
caxis([9.5 16.5]);
cb=colorbar;
cb.Ticks=10:16;
cb.TickLabels={'Isol. tethered conv.','Emb. tethered conv.','Isol. elevated conv.',...
    'Emb. elevated conv.','Stratiform','Strat. emb. conv.',...
    'Small'};
grid on
cb.Position=[0.8689     0.0500    0.0178    0.2098];
title('Stratiform/convective partitioning')

ax4.Position=[0.06 0.05 0.8 0.21];

linkaxes([ax1 ax2 ax3 ax4],'xy');

formatOut = 'yyyymmdd_HHMM';
set(gcf,'PaperPositionMode','auto')
print([figdir,'example',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
