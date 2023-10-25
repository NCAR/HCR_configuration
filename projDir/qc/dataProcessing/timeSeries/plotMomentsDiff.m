function plotMomentsDiff(data,momentsTime,momentsSpec,timeBeams,figdir,project,ylimUpper,flipYes,showPlot)
tFields=fields(momentsTime);
sFields=fields(momentsSpec);

commFields=intersect(tFields,sFields);

for ii=1:length(commFields)

f1 = figure('Position',[200 500 900 900],'DefaultAxesFontSize',12,'visible',showPlot);

colormap jet
s1=subplot(3,1,1);

hold on
surf(timeBeams,data.range./1000,momentsTime.(commFields{ii}),'edgecolor','none');
view(2);
ylabel('Range (km)');
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title([commFields{ii},' time domain']);

if flipYes
    set(gca, 'YDir','reverse');
end

s2=subplot(3,1,2);

hold on
surf(timeBeams,data.range./1000,momentsSpec.(commFields{ii}),'edgecolor','none');
view(2);
ylabel('Range (km)');
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title([commFields{ii},' spectral domain']);

if flipYes
    set(gca, 'YDir','reverse');
end

s3=subplot(3,1,3);

hold on
surf(timeBeams,data.range./1000,momentsSpec.(commFields{ii})-momentsTime.(commFields{ii}),'edgecolor','none');
view(2);
ylabel('Range (km)');
clim([-0.1 0.1]);
ylim([0 ylimUpper]);
xlim([timeBeams(1),timeBeams(end)]);
colorbar
grid on
title([commFields{ii},' spectral minus time domain']);

if flipYes
    set(gca, 'YDir','reverse');
end

if strcmp(commFields{ii},'dbz')
    s1.CLim=[-60,20];
    s2.CLim=[-60,20];
elseif strcmp(commFields{ii},'vel')
    s1.CLim=[-5,5];
    s2.CLim=[-5,5];
elseif strcmp(commFields{ii},'powerV')
    s1.CLim=[-110,-40];
    s2.CLim=[-110,-40];
elseif strcmp(commFields{ii},'width')
    s1.CLim=[0,4];
    s2.CLim=[0,4];
elseif strcmp(commFields{ii},'snr')
    s1.CLim=[-20,70];
    s2.CLim=[-20,70];
elseif strcmp(commFields{ii},'skew')
    s1.CLim=[-1,1];
    s2.CLim=[-1,1];
elseif strcmp(commFields{ii},'kurt')
    s1.CLim=[0,30];
    s2.CLim=[0,30];
end

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_momentsDiff_',(commFields{ii}),'_',datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0');
end
end