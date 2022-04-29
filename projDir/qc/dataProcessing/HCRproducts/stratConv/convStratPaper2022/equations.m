% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

showPlot='on';

blockTransition=1; % Remove data where antenna is in transition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

startTime=datetime(2019,9,27,12,39,0);
endTime=datetime(2019,9,27,12,44,0);

figdir=['/scr/sci/romatsch/other/convStratPaperHCR/'];

%% Get data

disp("Getting data ...");

fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ_MASKED=[];
data.VEL_MASKED=[];
data.TOPO=[];
data.TEMP=[];
data.MELTING_LAYER=[];
data.FLAG=[];

if blockTransition
    data.ANTFLAG=[];
end

dataVars=fieldnames(data);

% Load data
data=read_HCR(fileList,data,startTime,endTime);

% Check if all variables were found
for xInd2=1:length(dataVars)
    if ~isfield(data,dataVars{xInd2})
        dataVars{xInd2}=[];
    end
end

dataVars=dataVars(~cellfun('isempty',dataVars));

ylimUpper=(max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5;

% Take care of up pointing VEL
data.VEL_MASKED(:,data.elevation>0)=-data.VEL_MASKED(:,data.elevation>0);

% Remove data where antenna is in transition
if blockTransition
    data.DBZ_MASKED(:,data.ANTFLAG==5)=nan;
    data.VEL_MASKED(:,data.ANTFLAG==5)=nan;
end

%% Texture from reflectivity and velocity

disp('Calculating reflectivity texture ...');

pixRad=50; % Radius over which texture is calculated in pixels. Default is 50.
dbzBase=-10; % Reflectivity base value which is subtracted from DBZ.

dbzText=nan(size(data.DBZ_MASKED));

%DBZ(DBZ<dbzThresh)=nan;

% Pad data at start and end
dbzPadded=padarray(data.DBZ_MASKED,[0 pixRad],nan);

% Fill in areas with no data
dbzPadded=fillmissing(dbzPadded,'linear',2,'EndValues','nearest');

% Adjust reflectivity with base value
dbzPadded=dbzPadded-dbzBase;

% Loop through data points in time direction and pull out right window
%for ii=1:size(dbzPadded,2)-pixRad*2-1

%% Convective case
xInd=2590;
yInd=500;

dbzBlock=dbzPadded(:,xInd:xInd+pixRad*2);
timeBlock=data.time(xInd:xInd+pixRad*2);

timeInd=round(length(timeBlock)/2);

% Calculate and remove slope of reflectivity
% Calculate fit,'-b','LineWidth',1.5
x1=1:size(dbzBlock,2);
X=repmat(x1,size(dbzBlock,1),1);

sumX=sum(X,2,'omitnan');
sumY=sum(dbzBlock,2,'omitnan');
sumXY=sum((dbzBlock.*X),2,'omitnan');
sumX2=sum(X.^2,2,'omitnan');
sumY2=sum(dbzBlock.^2,2,'omitnan');

N=size(dbzBlock,2);

a=(sumY.*sumX2-sumX.*sumXY)./(N.*sumX2-sumX.^2);
b=(N.*sumXY-sumX.*sumY)./(N.*sumX2-sumX.^2);

newY=a+b.*X;

% Remove slope
dbzCorr=dbzBlock-newY+mean(dbzBlock,2,'omitnan');
dbzCorr(dbzCorr<1)=1;

% Calculate texture
tdbz=sqrt(std(dbzCorr.^2,[],2,'omitnan'));

%% Stratiform case
xInd2=650;
yInd2=110;

dbzBlock2=dbzPadded(:,xInd2:xInd2+pixRad*2);
timeBlock2=data.time(xInd2:xInd2+pixRad*2);

timeInd2=round(length(timeBlock2)/2);

% Calculate and remove slope of reflectivity
% Calculate fit,'-b','LineWidth',1.5
x1=1:size(dbzBlock2,2);
X=repmat(x1,size(dbzBlock2,1),1);

sumX=sum(X,2,'omitnan');
sumY=sum(dbzBlock2,2,'omitnan');
sumXY=sum((dbzBlock2.*X),2,'omitnan');
sumX2=sum(X.^2,2,'omitnan');
sumY2=sum(dbzBlock2.^2,2,'omitnan');

N=size(dbzBlock2,2);

a=(sumY.*sumX2-sumX.*sumXY)./(N.*sumX2-sumX.^2);
b=(N.*sumXY-sumX.*sumY)./(N.*sumX2-sumX.^2);

newY2=a+b.*X;

% Remove slope
dbzCorr2=dbzBlock2-newY2+mean(dbzBlock2,2,'omitnan');
dbzCorr2(dbzCorr2<1)=1;

% Calculate texture
tdbz2=sqrt(std(dbzCorr2.^2,[],2,'omitnan'));
    
    %% Plot

close all

wi=9;
hi=7;

fig1=figure('DefaultAxesFontSize',13,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
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
surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
view(2);
ylabel('Altitude (km)');
caxis([-35 25]);
ylim([0 ylimUpper]);
xlim([data.time(1),data.time(end)]);
cb1=colorbar;
plot(data.time,data.altitude./1000,'-k','LineWidth',2);
grid on
box on
text(datetime(2019,9,27,12,39,5),14.8,'(a) Reflectivity (dBZ)','FontSize',11,'FontWeight','bold');

plot([timeBlock(1),timeBlock(end)],[data.asl(yInd,timeInd)./1000,data.asl(yInd,timeInd)./1000],'-k','LineWidth',1.5);
scatter(timeBlock(timeInd),data.asl(yInd,timeInd)./1000,'filled','k')
text(timeBlock(timeInd)-seconds(15),data.asl(yInd,timeInd)./1000,'B','FontSize',13,'FontWeight','bold')

plot([timeBlock2(1),timeBlock2(end)],[data.asl(yInd2,timeInd2)./1000,data.asl(yInd2,timeInd2)./1000],'-k','LineWidth',1.5);
scatter(timeBlock2(timeInd2),data.asl(yInd2,timeInd2)./1000,'filled','k')
text(timeBlock2(timeInd)-seconds(15),data.asl(yInd2,timeInd2)./1000,'A','FontSize',13,'FontWeight','bold')

ax = gca;
ax.SortMethod='childorder';

s2=subplot(3,2,3);

hold on
l1=plot(timeBlock2,dbzBlock2(yInd2,:)+dbzBase,'-b','LineWidth',1.5);
l2=plot(timeBlock2,newY2(yInd2,:)+dbzBase,'-r','LineWidth',1.5);
xlim([timeBlock2(1),timeBlock2(end)]);
ylim([-15 5]);
grid on
box on
set(gca,'XTickLabel',[]);
ylabel('DBZ, DBZ_{fit}')
text(datetime(2019,9,27,12,40,5),6.5,'(b) A: DBZ and DBZ_{fit} (dBZ)','FontSize',11,'FontWeight','bold');
legend([l1,l2],{'DBZ','DBZ_{fit}'},'orientation','horizontal','Location','northwest')

s3=subplot(3,2,5);

hold on
plot(timeBlock2,dbzCorr2(yInd2,:)+dbzBase,'-b','LineWidth',1.5);
xlim([timeBlock2(1),timeBlock2(end)]);
ylim([-15 5]);
grid on
box on
ylabel('DBZ_{corr}')
text(timeBlock2(1)+seconds(0.2),-13.5,['TDBZ: ',num2str(tdbz2(yInd2),2),' dBZ, Scaled TDBZ: ',num2str(tdbz2(yInd2)/12,2)],'FontSize',12)
text(datetime(2019,9,27,12,40,5),6.5,'(d) A: DBZ_{corr} and DBZ_{adj} (dBZ)','FontSize',11,'FontWeight','bold');

yyaxis right
plot([0,0],[0,0]);
ax = gca;
ax.YColor='k';
ylim([-15-dbzBase 5-dbzBase]);
ylh=ylabel('DBZ_{adj}');
set(ylh,'rotation',-90,'VerticalAlignment','bottom')

s4=subplot(3,2,4);

hold on
l1=plot(timeBlock,dbzBlock(yInd,:)+dbzBase,'-b','LineWidth',1.5);
l2=plot(timeBlock,newY(yInd,:)+dbzBase,'-r','LineWidth',1.5);
xlim([timeBlock(1),timeBlock(end)]);
ylim([-5 15]);
set(gca,'XTickLabel',[]);
grid on
box on
ylabel('DBZ, DBZ_{fit}')
text(datetime(2019,9,27,12,43,19),16.5,'(c) B: DBZ and DBZ_{fit} (dBZ)','FontSize',11,'FontWeight','bold');
legend([l1,l2],{'DBZ','DBZ_{fit}'},'orientation','horizontal','Location','northwest')

s5=subplot(3,2,6);

hold on
plot(timeBlock,dbzCorr(yInd,:)+dbzBase,'-b','LineWidth',1.5);
xlim([timeBlock(1),timeBlock(end)]);
ylim([-5 15]);
grid on
box on
ylabel('DBZ_{corr}')
text(timeBlock(1)+seconds(0.2),-3.5,['TDBZ: ',num2str(tdbz(yInd),3),' dBZ, Scaled TDBZ: ',num2str(tdbz(yInd)/12,2)],'FontSize',12)
text(datetime(2019,9,27,12,43,19),16.5,'(e) B: DBZ_{corr} and DBZ_{adj} (dBZ)','FontSize',11,'FontWeight','bold');

yyaxis right
plot([0,0],[0,0]);
ax = gca;
ax.YColor='k';
ylim([-5-dbzBase 15-dbzBase]);
ylh=ylabel('DBZ_{adj}');
set(ylh,'rotation',-90,'VerticalAlignment','bottom')

s1.Position=[0.08 0.72 0.845 0.24];
cb1.Position=[0.93 0.72 0.0247 0.24];

s2.Position=[0.08 0.37 0.35 0.24];
s3.Position=[0.08 0.07 0.35 0.24];

s4.Position=[0.58 0.37 0.35 0.24];
s5.Position=[0.58 0.07 0.35 0.24];

set(gcf,'PaperPositionMode','auto')
print(fig1,[figdir,'equations.png'],'-dpng','-r0')
