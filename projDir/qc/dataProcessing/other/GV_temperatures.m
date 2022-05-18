% Compare HCR temperatures

clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir='/scr/sleet2/rsfdata/projects/otrec/GV/LRT/';

figdir=['/scr/sleet2/rsfdata/projects/otrec/GV/LRT/figs/'];

infiles=dir([indir,'RF*']);

for ii=1:length(infiles)
    file=[indir,infiles(ii).name];

    disp(file)

    baseTime=datetime(str2num(infiles(ii).name(6:9)),str2num(infiles(ii).name(10:11)),str2num(infiles(ii).name(12:13)));
    secondsIn=ncread(file,'Time');
    time=baseTime+seconds(secondsIn);
    temp=ncread(file,'ATX');

    close all
    
    f1=figure('DefaultAxesFontSize',14,'renderer','painters');
    set(f1,'Position',[200 500 1500 500]);
    hold on;
    plot(time,temp,'LineWidth',2);
    xlim([time(1),time(end)]);
    
    ylabel('Temperature [C]');

    grid on
    box on

    formatOut = 'yyyymmdd_HHMMSS';
    timestring=datestr(time(1),formatOut);
    timestring2=datestr(time(end),formatOut);
    title([timestring,' to ',timestring2],'interpreter','none');

    set(gcf,'PaperPositionMode','auto')
    print(f1, [figdir,'RF',num2str(ii),'_temperatures_',timestring,'_to_',timestring2],'-dpng','-r0');

end