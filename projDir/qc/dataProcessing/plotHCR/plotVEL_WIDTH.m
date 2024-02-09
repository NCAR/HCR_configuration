% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset, otrec
quality='qc1'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v1.2';

startTime=datetime(2021,6,11,19,36,0);
endTime=datetime(2021,6,11,19,37,30);

indir=HCRdir(project,quality,qcVersion,freqData);

ylims=[2,4.5];

saveFig=1;
if saveFig
    outname=[project,'_vel_width_',datestr(startTime,'yyyymmdd_HHMMSS')];
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

data.VEL_UNFOLDED=[];
data.WIDTH=[];

data=read_HCR(fileList,data,startTime,endTime);

data.VEL_UNFOLDED(:,data.elevation>=0)=-data.VEL_UNFOLDED(:,data.elevation>=0);
data.WIDTH(isnan(data.VEL_UNFOLDED))=nan;

%% Plot
disp('Plotting ...');

fig=figure('Position',[200 500 1200 1000],'DefaultAxesFontSize',14);

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);

surf(data.time,data.asl./1000,data.VEL_UNFOLDED,'EdgeColor','none');
view(2)
clim([-16 16]);
s1.Colormap=velCols;
colorbar

xlim([startTime,endTime]);
ylim([ylims]);

ylabel('Altitude (km)')

grid on
box on

title('Velocity (m s^{-1}) (positive is down)')

s2=nexttile(2);

surf(data.time,data.asl./1000,data.WIDTH,'EdgeColor','none');
view(2)
clim([0 4]);
s2.Colormap=turbo;
colorbar

xlim([startTime,endTime]);
ylim([ylims]);

ylabel('Altitude (km)')

grid on
box on

title('Spectrum width (m s^{-1})')

if saveFig
    set(gcf,'PaperPositionMode','auto')
    print(fig,[figdir,outname,'.png'],'-dpng','-r0')
end