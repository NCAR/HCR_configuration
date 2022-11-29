% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

startTime=datetime(2019,9,21,16,30,0);
endTime=datetime(2019,9,21,16,42,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=15;

saveFig=1;
if saveFig
    outname=['multScatt_',datestr(startTime,'yyyymmdd_HHMMSS')];
    %figdir=['/scr/sci/romatsch/HCR/examplePlots/'];
    figdir=[indir(1:end-5),'examplePlots/'];
    if ~exist(figdir, 'dir')
        mkdir(figdir)
    end
end

%% Get data

disp('Reading data ...');

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

data=[];

data.DBZ=[];
data.DBMVC=[];
data.DBMHX=[];
data.LDR=[];

data=read_HCR(fileList,data,startTime,endTime);

data.DBMHX(isnan(data.LDR))=-106;
data.DBMHX(isnan(data.DBZ))=nan;

data.DBMVC(isnan(data.DBZ))=nan;

%% Plot
close all
disp('Plotting ...');

plotIndNum=1000;
plotSpace=round(length(data.time)/plotIndNum);
plotInds=1:plotSpace:length(data.time);

fig=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',14);

colormap('jet');

s1=subplot(3,1,1);

surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.DBMVC(:,plotInds),'EdgeColor','none');
view(2)
caxis([-105,-70]);
colorbar

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

title('DBMVC (dB)')

s2=subplot(3,1,2);
colmapHX=jet(63);
colmapHX=cat(1,[0.7,0,0],colmapHX);

surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.DBMHX(:,plotInds),'EdgeColor','none');
view(2)
caxis([-106,-80]);
colormap(s2,colmapHX);
colorbar

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

title('DBMHX (dB)')

s3=subplot(3,1,3);

surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.LDR(:,plotInds),'EdgeColor','none');
view(2)
caxis([-25,0]);
colorbar

xlim([startTime,endTime]);
ylim([0 ylimUpper]);

ylabel('Altitude (km)')

grid on
box on

title('LDR (dB)')

linkaxes([s1 s2 s3],'xy')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end