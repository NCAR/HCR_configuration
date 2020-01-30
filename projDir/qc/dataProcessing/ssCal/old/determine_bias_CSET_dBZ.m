% Ocean scan calibration for HCR data

clear all;
close all;

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made
makeSingleFigs=0;

addpath('/h/eol/romatsch/git/private/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/git/private/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/git/private/utils/');

directories.figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/csetRaw/';
%directories.figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/csetQC/';

directories.dataDir='/scr/eldora2/rsfdata/cset/hcr/cfradial/moments/10hz/'; %raw data
%directories.dataDir='/scr/eldora2/rsfdata/cset/hcr/qc/cfradial/moments/10hz/'; %qc data

directories.sondedir='/h/eol/romatsch/data/hcrCalib/soundings/cset_drops/';
%directories.modeldir='/h/eol/romatsch/data/reanalysis/ecmwf/socrates/soundIn/';

infile='/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/cal_CSET.txt';

projectName='CSET';

[outTable PLT]=f_processOceanScans(directories,infile,makeSingleFigs,projectName);


%% Plot bias vs wind speed
close all
figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(outTable.SfcWindSpd,outTable.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
set(gca,'XLim',[0 15]);
set(gca,'YLim',[0,7]);
xlabel('Surface wind speed [m/s]');
ylabel('Liebe bias [dB]');
title([projectName ' Liebe bias vs surface wind speed']);
plabel = [directories.figdir projectName '_biasVSwindspd_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot bias std vs wind speed
close all
figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(outTable.SfcWindSpd,outTable.Std, 'o', 'MarkerFaceColor', 'b');
text(outTable.SfcWindSpd,outTable.Std+0.05,num2str(outTable.N),'HorizontalAlignment','center');
set(gca,'XLim',[0 15]);
set(gca,'YLim',[0,1.5]);
xlabel('Surface wind speed [m/s]');
ylabel('Liebe bias standard deviation [dB]');
title([projectName ' Liebe bias standard deviation vs surface wind speed']);
plabel = [directories.figdir projectName '_biasStcVSwindspd_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot bias vs wind angle

hdg=outTable.MeanHdg;
hdg(hdg<0)=hdg(hdg<0)+360;

normDeg = mod(hdg-outTable.SfcWindSpd,360);
hdgMinusDir=min(360-normDeg, normDeg);

close all
f2 = figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(hdgMinusDir,outTable.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
set(gca,'XLim',[0 180]);
set(gca,'YLim',[0,7]);
xlabel('abs(heading - surface wind dir) [deg]');
ylabel('Liebe bias [dB]');
title([projectName ' Liebe bias vs heading minus wind direction']);
plabel = [directories.figdir projectName '_biasVSwinddir_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

