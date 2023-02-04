% Analyze HCR clouds

clear all;
close all;

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.1';
freqData='10hz'; % 10hz, 100hz, or 2hz

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc3/cfradial/v3.1_full/';

indir=HCRdir(project,quality,qcVersion,freqData);

%% Load data1

disp('Loading data ...');

startTime=datetime(2019,8,11,13,53,0);
endTime=datetime(2019,8,11,14,26,0);

data=[];

data.MELTING_LAYER=[];
data.ICING_LEVEL=[];
data.VEL_MASKED=[];
data.DBZ_MASKED=[];
data.ECHO_TYPE_2D=[];
data.PID=[];

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

elevenInds=find(data.MELTING_LAYER==11);
twelveInds=find(data.MELTING_LAYER==12);
thirteenInds=find(data.MELTING_LAYER==13);
fourteenInds=find(data.MELTING_LAYER==14);

twentyoneInds=find(data.MELTING_LAYER==21);
twentytwoInds=find(data.MELTING_LAYER==22);
twentythreeInds=find(data.MELTING_LAYER==23);
twentyfourInds=find(data.MELTING_LAYER==24);

timeMat=repmat(data.time,size(data.MELTING_LAYER,1),1);

if etime(datevec(endTime),datevec(startTime))<=900
    newInds=1:1:length(data.time);
elseif etime(datevec(endTime),datevec(startTime))<=3600
    newInds=1:10:length(data.time);
else
    newInds=1:100:length(data.time);
end

% Resample for plotting
newFindMelt=data.MELTING_LAYER(:,newInds);
newASL=data.asl(:,newInds);
newVEL=data.VEL_MASKED(:,newInds);
newDBZ=data.DBZ_MASKED(:,newInds);
newPID=data.PID(:,newInds);
newTime=data.time(newInds);

% Echo type
stratConvPlot=data.ECHO_TYPE_2D(:,newInds);
stratConvPlot(stratConvPlot==14)=1;
stratConvPlot(stratConvPlot==16)=2;
stratConvPlot(stratConvPlot==18)=3;
stratConvPlot(stratConvPlot==25)=4;
stratConvPlot(stratConvPlot==30)=5;
stratConvPlot(stratConvPlot==32)=6;
stratConvPlot(stratConvPlot==34)=7;
stratConvPlot(stratConvPlot==36)=8;
stratConvPlot(stratConvPlot==38)=9;

colmapSC=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,0;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    0.7,0,0];

cscale_hcr=[255,0,0; 255,204,204; 249,163,25; 255,240,60; 136,34,185; 255,0,255; 17,170,51; 0,0,255; 0,255,255; 0,0,0; 150,150,150];

cscale_hcr=cscale_hcr./255;

units_str_hcr={'Rain','SC Rain','Drizzle','SC Drizzle','Cloud Liquid','SC Cloud Liq.',...
    'Melting','Large Frozen','Small Frozen','Precip','Cloud'};
%% Plot
close all

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

ax1=subplot(4,1,1);
hold on;
surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
ax1.Colormap=dbz_default;
view(2);
caxis([-60 25]);
scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,5,'g','filled','filled','MarkerEdgeColor','g');

ax1.SortMethod = 'childorder';
ylim([0 14]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Reflectivity (dBZ) and melting layer (green)')
grid on
box on
cb1=colorbar;
set(gca,'XTickLabel',[]);

ax2=subplot(4,1,2);
hold on;
surf(newTime,newASL./1000,newVEL,'edgecolor','none');
caxis([-6 6]);
cb2=colorbar;
ax2.Colormap=velCols;
cb2.Ticks=-6:2:6;
view(2);

ylim([0 14]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Doppler velocity (m s^{-1})')
grid on
box on
set(gca,'XTickLabel',[]);

ax3=subplot(4,1,3);
hold on;
surf(newTime,newASL./1000,stratConvPlot,'edgecolor','none');
view(2);

ylim([0 14]);
ylabel('Altitude (km)');
xlim([data.time(1),data.time(end)]);
title('Convective/stratiform echo type')
grid on
box on
ax3.Colormap=colmapSC;
caxis([0.5 9.5]);
cb3=colorbar;
cb3.Ticks=1:9;
cb3.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
    'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
set(gca,'XTickLabel',[]);

ax4=subplot(4,1,4);
surf(newTime,newASL./1000,newPID,'edgecolor','none');
view(2);
ylim([0 14]);
xlim([data.time(1),data.time(end)]);
caxis([0.5 11.5]);
colormap(ax4,cscale_hcr);
cb4=colorbar;
cb4.Ticks=1:11;
cb4.TickLabels=units_str_hcr;
ylabel('Altitude (km)');
title('Particle identifications')
grid on
box on

ax1.Position=[0.05,0.765,0.8,0.195];
ax2.Position=[0.05,0.53,0.8,0.195];
ax3.Position=[0.05,0.295,0.8,0.195];
ax4.Position=[0.05,0.06,0.8,0.195];

print([figdir,'HCRvariables.png'],'-dpng','-r0');

