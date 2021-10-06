% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

minPixNumUW=5;

HCRrangePix=10;
HCRtimePix=20;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';

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

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/paperFigs/'];

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
units_str_hcr={'Rain','SC rain','Drizzle','SC drizzle','Cloud liquid','SC cloud liquid','Mixed phase','Large frozen','Small frozen'};

cscale_hcr_2=[1 0 0;0 1 0;0 0 1];
units_str_hcr_2={'Liquid','Mixed','Frozen'};

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

data.DBZ(data.FLAG>1)=nan;

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
hcrLiqIce(data.PID<=6)=1;
hcrLiqIce(data.PID==7)=2;
hcrLiqIce(data.PID>=8)=3;

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

%% Plot 1


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

wi=8;
hi=6;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(3,1,1);

colormap jet

hold on
surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
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

s2=subplot(3,1,2);

hold on
surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
view(2);
colormap(s2,cscale_hcr);
cb2=colorbar;
cb2.Ticks=1:9;
cb2.TickLabels=units_str_hcr;
ylabel('Altitude (km)');
%title(['HCR particle ID']);
text(startTime+seconds(60),ylims(2)-0.11,'(b) PID',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);

l1=plot(data.time,data.altitude./1000,'-k','linewidth',2);
legend(l1,'Aircraft altitude');

caxis([.5 9.5]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);
s2.XTickLabel=[];

grid on
box on

% Dummy for colorbar
s4=subplot(4,1,4);

colormap(s4,colmapL(2:end,:));
cb4=colorbar;
cb4.Title.String='Liq. frac.       ';
cb4.Ticks=[0 0.5 1];
s4.Position=[5 0.451 0.79 0.13];

s3=subplot(3,1,3);

hold on
surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
view(2);
colormap(s3,cscale_hcr_2);
cb3=colorbar;
cb3.Ticks=6:8;
cb3.TickLabels=units_str_hcr_2;
cb3.Title.String=' PID';
ylabel('Altitude (km)');
%title(['HCR particle ID']);
text(startTime+seconds(60),ylims(2)-0.11,'(c) Simplified PID and UWILD liquid fraction',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);

scatter(ptime,ttSync.Var1./1000,20,col1DL,'filled');
set(gca,'clim',[0,1]);

caxis([5.5 8.5]);
ylim(ylims);
xlim([data.time(1),data.time(end)]);

grid on
box on

s1.Position=[0.065 0.7 0.76 0.28];
cb1.Position=[0.84 0.7 0.023 0.28];
cb1.FontSize=9;

s2.Position=[0.065 0.39 0.76 0.28];
cb2.Position=[0.84 0.39 0.023 0.28];
cb2.FontSize=9;

s3.Position=[0.065 0.08 0.76 0.28];
cb3.Position=[0.91 0.08 0.023 0.25];
cb3.FontSize=9;

cb4.Position=[0.87 0.08 0.023 0.25];
cb4.FontSize=9;

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'pidUW.png'],'-dpng','-r0')
