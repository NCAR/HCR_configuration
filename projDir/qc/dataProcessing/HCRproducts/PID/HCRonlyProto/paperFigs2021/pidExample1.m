% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

ylimits=[0 2];

indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';

figdir=['/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/pidPlots/paperFigs/'];

startTime=datetime(2018,1,19,4,25,0);
endTime=datetime(2018,1,19,4,30,0);

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

%% Load data

disp('Loading data');

data=[];

%HCR data
data.DBZ=[];
data.VEL_CORR=[];
data.WIDTH=[];
data.LDR=[];
data.TEMP=[];
data.MELTING_LAYER=[];
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

% Mask with FLAG
data.DBZ(data.FLAG>1)=nan;
data.VEL_CORR(data.FLAG>1)=nan;
data.WIDTH(data.FLAG>1)=nan;
data.LDR(data.FLAG>1)=nan;
%data.TEMP(data.FLAG>1)=nan;
%data.MELTING_LAYER(data.FLAG>1)=nan;

%% Scales and units
cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];

units_str_hcr={'Rain','SC rain','Drizzle','SC drizzle','Cloud liquid','SC cloud liquid','Mixed phase','Large frozen','Small frozen'};

%% Plot PIDs

timeMat=repmat(data.time,size(data.TEMP,1),1);

disp('Plotting PID');

close all

wi=8;
hi=9.5;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(7,1,1);
surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-40 20]);
colormap(s1,jet);
cb1=colorbar;
cb1.Ticks=-30:10:10;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.2,'(a) DBZ (dBZ)',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on
s1.XTickLabel=[];
s1.YTick=0:0.5:1.5;

s2=subplot(7,1,2);
surf(data.time,data.asl./1000,data.VEL_CORR,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-5 5]);
colormap(s2,jet);
cb2=colorbar;
cb2.Ticks=-4:2:4;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.25,'(b) VEL (m s^{-1})',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
s2.SortMethod = 'childorder';
grid on
box on
s2.XTickLabel=[];
s2.YTick=0:0.5:1.5;

s3=subplot(7,1,3);
surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-30 -5]);
colormap(s3,jet);
cb3=colorbar;
cb3.Ticks=-25:5:-10;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.2,'(c) LDR (dB)',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on
s3.XTickLabel=[];
s3.YTick=0:0.5:1.5;

s4=subplot(7,1,4);
surf(data.time,data.asl./1000,data.WIDTH,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([0 2]);
colormap(s4,jet);
cb4=colorbar;
cb4.Ticks=0.25:0.25:1.75;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.25,'(d) WIDTH (m s^{-1})',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
s4.SortMethod = 'childorder';
grid on
box on
s4.XTickLabel=[];
s4.YTick=0:0.5:1.5;

s5=subplot(7,1,5);
jetIn=jet;
jetTemp=cat(1,jetIn(1:size(jetIn,1)/2,:),repmat([0 0 0],3,1),...
    jetIn(size(jetIn,1)/2+1:end,:));
surf(data.time,data.asl./1000,data.TEMP,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([-10 10]);
colormap(s5,jetTemp);
cb5=colorbar;
cb5.Ticks=-8:2:8;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.2,['(e) TEMP ',char(176),'C'],...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
s5.SortMethod = 'childorder';
grid on
box on
s5.XTickLabel=[];
s5.YTick=0:0.5:1.5;

s6=subplot(7,1,6);
plotMelt=data.MELTING_LAYER;
plotMelt(~isnan(plotMelt) & plotMelt<20)=10;
plotMelt(~isnan(plotMelt) & plotMelt>=20)=20;
surf(data.time,data.asl./1000,plotMelt,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
%caxis([0 2]);
colormap(s6,[1 0 1;1 1 0]);
cb6=colorbar;
cb6.Ticks=[12.5 17.5];
cb6.TickLabels={'Below ICING','Above ICING'};
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.2,'(f) ICING',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
s6.SortMethod = 'childorder';
grid on
box on
s6.XTickLabel=[];
s6.YTick=0:0.5:1.5;

s7=subplot(7,1,7);
surf(data.time,data.asl./1000,data.PID,'edgecolor','none');
view(2);
ylim(ylimits);
xlim([data.time(1),data.time(end)]);
caxis([.5 9.5]);
colormap(s7,cscale_hcr);
cb7=colorbar;
cb7.Ticks=1:9;
cb7.TickLabels=units_str_hcr;
ylabel('Altitude (km)');
text(startTime+seconds(5),ylimits(2)-0.2,'(g) PID',...
    'fontsize',11,'fontweight','bold','BackgroundColor','w','Margin',0.5);
grid on
box on
s7.YTick=0:0.5:1.5;

s1.Position=[0.057 0.862 0.79 0.13];
cb1.Position=[0.852 0.862 0.023 0.13];
cb1.FontSize=9;

s2.Position=[0.057 0.725 0.79 0.13];
cb2.Position=[0.852 0.725 0.023 0.13];
cb2.FontSize=9;

s3.Position=[0.057 0.588 0.79 0.13];
cb3.Position=[0.852 0.588 0.023 0.13];
cb3.FontSize=9;

s4.Position=[0.057 0.451 0.79 0.13];
cb4.Position=[0.852 0.451 0.023 0.13];
cb4.FontSize=9;

s5.Position=[0.057 0.314 0.79 0.13];
cb5.Position=[0.852 0.314 0.023 0.13];
cb5.FontSize=9;

s6.Position=[0.057 0.177 0.79 0.13];
cb6.Position=[0.852 0.177 0.023 0.13];
cb6.FontSize=9;

s7.Position=[0.057 0.04 0.79 0.13];
cb7.Position=[0.852 0.04 0.023 0.13];
cb7.FontSize=9;

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'pidExample1.png'],'-dpng','-r0')
