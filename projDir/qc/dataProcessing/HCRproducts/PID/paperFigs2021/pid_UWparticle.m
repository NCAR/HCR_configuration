% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%
project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='combined'; % 10hz, 100hz, 2hz, or combined
%whichModel='era5';

largeUW=1; % Set to 1 when we want to use only larges particles up to minPixNumUW
minPixNumUW=20;

coldOnly=0; % Set to 1 when only cold region is desired

HCRrangePix=5;
HCRtimePix=4;
minPixNumHCR=14;

plotOn=0;
showPlot='off';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

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
figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc3/cfradial/hcr_hsrl_merge/v3.0_full/pidPlotsComb/paperFigs/'];

%cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
cscale_hcr=[255,0,0; 255,204,204; 249,163,25; 255,240,60; 136,34,185; 255,0,255; 17,170,51; 0,0,255; 0,255,255; 0,0,0; 150,150,150];
cscale_hcr=cscale_hcr./255;
units_str_hcr={'Rain','SC Rain','Drizzle','SC Drizzle','Cloud Liquid','SC Cloud Liq.',...
    'Melting','Large Frozen','Small Frozen','Precip','Cloud'};

cscale_hcr_2=[1 0 0;0 1 0;0 0 1];
units_str_hcr_2={'Frozen','Melting','Liquid'};

varNames={'numLiqHCR','numIceHCR','numAllHCR','pidHCR','numLiqLargestP','numIceLargestP','numAllLargestP','sizeLargestP'};

aa=5;

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

startTime=datetime(2018,1,26,2,14,30);
endTime=datetime(2018,1,26,2,38,0);

%% Sub times for particles
ptimeInds=find(ptime1>=startTime & ptime1<=endTime);

countAll=countAll1(:,ptimeInds);

ptime=ptime1(ptimeInds);
countLiq=countLiq1(:,ptimeInds);
countIce=countIce1(:,ptimeInds);

sumAll=sum(countAll,1);

liqFrac=sum(countLiq,1)./sumAll;

%% Get HCR data

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.HCR_DBZ = [];
data.PID=[];
data.HSRL_Aerosol_Backscatter_Coefficient=[];

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

%% Find largest
countAllFlip=flipud(countAll);
cumSumAll=cumsum(countAllFlip);
cumSumAllBack=flipud(cumSumAll);

indMat=repmat((1:size(countAll,1))',1,size(countAll,2));
indMat(countAll==0)=nan;
indMat(cumSumAllBack<minPixNumUW)=nan;

rowsGood=max(indMat,[],1,'omitnan');
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
hcrLiqIce(data.PID<=6.1)=1;
hcrLiqIce(data.PID>=6.9 & data.PID<=7.1)=2;
hcrLiqIce(data.PID>=7.9 & data.PID<9.8)=3;

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

        if allNum>minPixNumHCR
            hcrPID=data.PID(18:18+HCRrangePix,hcrIndCols);
            pidOut=mode(reshape(hcrPID,1,[]));
            addFrac=[liqNum/allNum,liqFrac(goodIndsP(jj)),liqNum,iceNum,allNum,pidOut];
            liqFrac_HCR_P((goodIndsP(jj)),:)=addFrac;
        end
    end
end

outTable=timetable(ptime,liqFrac_HCR_P(:,3),liqFrac_HCR_P(:,4),liqFrac_HCR_P(:,5),liqFrac_HCR_P(:,6),...
    numLiqLargest',numIceLargest',numAllLargest',sizeLargest','VariableNames',varNames);

%% Plot 1
pidSimp=nan(size(hcrLiqIce));
pidSimp(hcrLiqIce==1)=3;
pidSimp(hcrLiqIce==3)=1;
pidSimp(hcrLiqIce==2)=2;

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

ylims=[0 1.6];

wi=10;
hi=9.3;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(4,1,1);

colormap jet

hold on
surf(data.time,data.asl./1000,data.HCR_DBZ,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);
s1.XTickLabel=[];
cb1=colorbar;
grid on
box on
text(startTime+seconds(60),ylims(2)-0.11,'(a) DBZ (dBZ)',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
%title('Reflectivity (dBZ)')
l1=plot(data.time,data.altitude./1000,'-k','linewidth',2);
legend(l1,'Aircraft altitude');

s2=subplot(4,1,2);

colormap jet

data.HSRL_Aerosol_Backscatter_Coefficient(isnan(data.HCR_DBZ))=nan;
data.HSRL_Aerosol_Backscatter_Coefficient(data.HSRL_Aerosol_Backscatter_Coefficient<9.9e-9)=nan;

hold on
surf(data.time,data.asl./1000,log10(data.HSRL_Aerosol_Backscatter_Coefficient),'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-9 0]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);
s2.XTickLabel=[];
cb2=colorbar;
grid on
box on
text(startTime+seconds(60),ylims(2)-0.11,'(b) log10(HSRL BACKSCAT) (m^{-1} sr^{-1})',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
%title('Reflectivity (dBZ)')
l1=plot(data.time,data.altitude./1000,'-k','linewidth',2);
legend(l1,'Aircraft altitude');

s3=subplot(4,1,3);

hold on
surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
view(2);
colormap(s3,cscale_hcr);
cb3=colorbar;
cb3.Ticks=1:11;
cb3.TickLabels=units_str_hcr;
ylabel('Altitude (km)');
%title(['HCR particle ID']);
text(startTime+seconds(60),ylims(2)-0.11,'(c) PID',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);

l1=plot(data.time,data.altitude./1000,'-k','linewidth',2);
legend(l1,'Aircraft altitude');

caxis([.5 11.5]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);
s3.XTickLabel=[];

grid on
box on

% Dummy for colorbar
s5=subplot(5,1,5);

colormap(s5,colmapL(2:end,:));
cb5=colorbar;
cb5.Title.String='Liq. frac.       ';
cb5.Ticks=[0 0.5 1];
s5.Position=[5 0.451 0.79 0.13];

s4=subplot(4,1,4);

hold on
surf(data.time,data.asl./1000,pidSimp,'edgecolor','none');
view(2);
colormap(s4,flipud(cscale_hcr_2));
cb4=colorbar;
cb4.Ticks=1:3;
cb4.TickLabels=units_str_hcr_2;
cb4.Title.String=' PID';
ylabel('Altitude (km)');
%title(['HCR particle ID']);
text(startTime+seconds(60),ylims(2)-0.11,'(d) Simplified PID and UWILD liquid fraction',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);

scatter(ptime,ttSync.Var1./1000,20,col1DL,'filled');
set(gca,'clim',[0,1]);

caxis([0.5 3.5]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);

grid on
box on

s1.Position=[0.065 0.767 0.775 0.23];
cb1.Position=[0.85 0.767 0.023 0.23];
cb1.FontSize=9;

s2.Position=[0.065 0.527 0.775 0.23];
cb2.Position=[0.85 0.527 0.023 0.23];
cb2.FontSize=9;

s3.Position=[0.065 0.287 0.775 0.23];
cb3.Position=[0.85 0.287 0.023 0.23];
cb3.FontSize=9;

s4.Position=[0.065 0.047 0.775 0.23];
cb4.Position=[0.92 0.047 0.023 0.21];
cb4.FontSize=9;

cb5.Position=[0.88 0.047 0.023 0.21];
cb5.FontSize=9;

set(gcf,'PaperPositionMode','auto')
%print(fig1,[figdir,'pidUW.png'],'-dpng','-r0')
print(fig1,[figdir,'pidUW.tif'],'-dtiffn','-r0')