% Ocean scan calibration for HCR data

clear all;
close all;

% If 1, plots for individual calibration events will be made, if 0, only
% total plots will be made
makeSingleFigs=1;

addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/functions/');
addpath('/h/eol/romatsch/gitPriv/process_HCR/oceanScans/colormaps/');
addpath('/h/eol/romatsch/gitPriv/utils/');

directories.figdir='/h/eol/romatsch/hcrCalib/oceanScans/figs/socratesRawWithECMWF/';

%directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/qc/cfradial/moments/10hz/'; %Final
directories.dataDir='/scr/rain1/rsfdata/projects/socrates/hcr/cfradial/moments/10hz/'; % raw data

directories.sondedir='/h/eol/romatsch/data/hcrCalib/soundings/socrates_drops/';
directories.modeldir='/scr/sci/romatsch/data/reanalysis/ecmwf/socrates/soundIn/';
%directories.wavedir='/scr/sci/romatsch/data/socrates_waves/netcdf/';

% Fast scans
infileF='/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/cal_SOCRATES_fast.txt';
projectNameF='fastSOCRATES';

[outTableFast PLTFast]=f_processOceanScans_dBZ(directories,infileF,makeSingleFigs,projectNameF);

% Slow scans
infileS='/h/eol/romatsch/hcrCalib/oceanScans/biasInFiles/cal_SOCRATES_slow.txt';
projectNameS='slowSOCRATES';

[outTableSlow PLTSlow]=f_processOceanScans_dBZ(directories,infileS,makeSingleFigs,projectNameS);

%% Merge fast and slow
PLTAll=cat(1,PLTFast,PLTSlow);

% Read file with calib events
caseListF = readtable(infileF,'Delimiter','space');
%convert to cell so each case has one cell
casesInF=table2array(caseListF(:,1));
numCasesF=unique(casesInF);
uniqueCasesF=cell(size(numCasesF,1),1);

for ii=1:size(numCasesF,1)
    caseIndF=find(casesInF==ii);
    uniqueCasesF{ii}=caseListF(caseIndF,:);
end

% Read file with calib events
caseListS = readtable(infileS,'Delimiter','space');
%convert to cell so each case has one cell
casesInS=table2array(caseListS(:,1));
numCasesS=unique(casesInS);
uniqueCasesS=cell(size(numCasesS,1),1);

for ii=1:size(numCasesS,1)
    caseIndS=find(casesInS==ii);
    uniqueCasesS{ii}=caseListS(caseIndS,:);
end

uniqueCasesA=cat(1,uniqueCasesF,uniqueCasesS);

reflXlims=[0 25];
reflYlims=[15,55];
sigXlims=[0 25];
sigYlims=[-20 15];

projectNameA='allSOCRATES';

plot_all_combined(PLTAll,sigXlims,sigYlims,reflXlims,reflYlims,projectNameA,directories,uniqueCasesA)


%% Plot bias vs wind speed
close all
figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(outTableFast.SfcWindSpd,outTableFast.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
plot(outTableSlow.SfcWindSpd,outTableSlow.BiasLiebe, 'oc', 'MarkerFaceColor', 'c');
set(gca,'XLim',[0 20]);
set(gca,'YLim',[0,4]);
xlabel('Surface wind speed [m/s]');
ylabel('Liebe bias [dB]');
title([projectNameA ' Liebe bias vs surface wind speed']);
legend('Fast','Slow');
plabel = [directories.figdir projectNameA '_biasVSwindspd_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot bias std vs wind speed
close all
figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(outTableFast.SfcWindSpd,outTableFast.Std, 'o', 'MarkerFaceColor', 'b');
plot(outTableSlow.SfcWindSpd,outTableSlow.Std, 'oc', 'MarkerFaceColor', 'c');
text(outTableFast.SfcWindSpd,outTableFast.Std-0.05,num2str(outTableFast.N),'HorizontalAlignment','center','color','b');
text(outTableSlow.SfcWindSpd,outTableSlow.Std+0.05,num2str(outTableSlow.N),'HorizontalAlignment','center','color','c');
set(gca,'XLim',[0 20]);
set(gca,'YLim',[0,1.5]);
xlabel('Surface wind speed [m/s]');
ylabel('Liebe bias standard deviation [dB]');
title([projectNameA ' Liebe bias standard deviation vs surface wind speed']);
legend('Fast','Slow');
plabel = [directories.figdir projectNameA '_biasStdVSwindspd_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot bias vs wind angle

hdgF=outTableFast.MeanHdg;
hdgF(hdgF<0)=hdgF(hdgF<0)+360;

normDegF = mod(hdgF-outTableFast.SfcWindSpd,360);
hdgMinusDirF=min(360-normDegF, normDegF);

hdgS=outTableSlow.MeanHdg;
hdgS(hdgS<0)=hdgS(hdgS<0)+360;

normDegS = mod(hdgS-outTableSlow.SfcWindSpd,360);
hdgMinusDirS=min(360-normDegS, normDegS);

close all
f2 = figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(hdgMinusDirF,outTableFast.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
plot(hdgMinusDirS,outTableSlow.BiasLiebe, 'oc', 'MarkerFaceColor', 'c');
set(gca,'XLim',[0 180]);
set(gca,'YLim',[0,4]);
xlabel('abs(heading - surface wind dir) [deg]');
ylabel('Liebe bias [dB]');
title([projectNameA ' Liebe bias vs heading minus wind direction']);
legend('Fast','Slow');
plabel = [directories.figdir projectNameA '_biasVSwinddir_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot std bias vs wind angle

hdgF=outTableFast.MeanHdg;
hdgF(hdgF<0)=hdgF(hdgF<0)+360;

normDegF = mod(hdgF-outTableFast.SfcWindSpd,360);
hdgMinusDirF=min(360-normDegF, normDegF);

hdgS=outTableSlow.MeanHdg;
hdgS(hdgS<0)=hdgS(hdgS<0)+360;

normDegS = mod(hdgS-outTableSlow.SfcWindSpd,360);
hdgMinusDirS=min(360-normDegS, normDegS);

close all
f2 = figure;
set(gcf,'Position',[200 500 800 600]);
hold on
plot(hdgMinusDirF,outTableFast.Std, 'o', 'MarkerFaceColor', 'b');
plot(hdgMinusDirS,outTableSlow.Std, 'oc', 'MarkerFaceColor', 'c');
set(gca,'XLim',[0 180]);
set(gca,'YLim',[0,1.5]);
xlabel('abs(heading - surface wind dir) [deg]');
ylabel('Liebe bias standard deviation [dB]');
title([projectNameA ' Liebe bias standard deviation vs heading minus wind direction']);
legend('Fast','Slow');
plabel = [directories.figdir projectNameA '_biasStdVSwinddir_all'];
set(gcf,'PaperPositionMode','auto')
print(plabel,'-dpng','-r0')

%% Plot bias std vs wave height
if isfield(directories,'wavedir')
    
    close all
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(outTableFast.SigWaveHt,outTableFast.Std, 'o', 'MarkerFaceColor', 'b');
    plot(outTableSlow.SigWaveHt,outTableSlow.Std, 'oc', 'MarkerFaceColor', 'c');
    text(outTableFast.SigWaveHt,outTableFast.Std-0.05,num2str(outTableFast.N),'HorizontalAlignment','center','color','b');
    text(outTableSlow.SigWaveHt,outTableSlow.Std+0.05,num2str(outTableSlow.N),'HorizontalAlignment','center','color','c');
    set(gca,'XLim',[0 8]);
    set(gca,'YLim',[0,1.5]);
    xlabel('Wave height [m]');
    ylabel('Liebe bias standard deviation [dB]');
    title([projectNameA ' Liebe bias standard deviation vs sig wave height']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasStdVSwaveHeight_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias vs wave height
    close all
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(outTableFast.SigWaveHt,outTableFast.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
    plot(outTableSlow.SigWaveHt,outTableSlow.BiasLiebe, 'oc', 'MarkerFaceColor', 'c');
    text(outTableFast.SigWaveHt,outTableFast.BiasLiebe-0.12,num2str(outTableFast.N),'HorizontalAlignment','center','color','b');
    text(outTableSlow.SigWaveHt,outTableSlow.BiasLiebe+0.12,num2str(outTableSlow.N),'HorizontalAlignment','center','color','c');
    set(gca,'XLim',[0 8]);
    set(gca,'YLim',[0,5]);
    xlabel('Wave height [m]');
    ylabel('Liebe bias [dB]');
    title([projectNameA ' Liebe bias vs sig wave height']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasVSwaveHeight_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias std vs swell wave height
    close all
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(outTableFast.SwellWaveHt,outTableFast.Std, 'o', 'MarkerFaceColor', 'b');
    plot(outTableSlow.SwellWaveHt,outTableSlow.Std, 'oc', 'MarkerFaceColor', 'c');
    text(outTableFast.SwellWaveHt,outTableFast.Std-0.05,num2str(outTableFast.N),'HorizontalAlignment','center','color','b');
    text(outTableSlow.SwellWaveHt,outTableSlow.Std+0.05,num2str(outTableSlow.N),'HorizontalAlignment','center','color','c');
    set(gca,'XLim',[0 5]);
    set(gca,'YLim',[0,1.5]);
    xlabel('Wave height [m]');
    ylabel('Liebe bias standard deviation [dB]');
    title([projectNameA ' Liebe bias standard deviation vs swell wave height']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasStdVSswellWaveHeight_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias vs wave height
    close all
    figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(outTableFast.SwellWaveHt,outTableFast.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
    plot(outTableSlow.SwellWaveHt,outTableSlow.BiasLiebe, 'oc', 'MarkerFaceColor', 'c');
    text(outTableFast.SwellWaveHt,outTableFast.BiasLiebe-0.12,num2str(outTableFast.N),'HorizontalAlignment','center','color','b');
    text(outTableSlow.SwellWaveHt,outTableSlow.BiasLiebe+0.12,num2str(outTableSlow.N),'HorizontalAlignment','center','color','c');
    set(gca,'XLim',[0 5]);
    set(gca,'YLim',[0,5]);
    xlabel('Wave height [m]');
    ylabel('Liebe bias [dB]');
    title([projectNameA ' Liebe bias vs swell wave height']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasVSswellWaveHeight_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot bias vs wave angle
    
    WnormDegF = mod(hdgF-outTableFast.SwellWaveDir,360);
    hdgMinusWDirF=min(360-WnormDegF, WnormDegF);
    
    WnormDegS = mod(hdgS-outTableSlow.SwellWaveDir,360);
    hdgMinusWDirS=min(360-WnormDegS, WnormDegS);
    
    close all
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(hdgMinusWDirF,outTableFast.BiasLiebe, 'o', 'MarkerFaceColor', 'b');
    plot(hdgMinusWDirS,outTableSlow.BiasLiebe, 'oc', 'MarkerFaceColor', 'c');
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0,4]);
    xlabel('abs(heading - swell wave dir) [deg]');
    ylabel('Liebe bias [dB]');
    title([projectNameA ' Liebe bias vs heading minus swell wave direction']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasVSwsbrdir_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
    %% Plot std bias vs wave angle
    
    WnormDegF = mod(hdgF-outTableFast.SwellWaveDir,360);
    hdgMinusWDirF=min(360-WnormDegF, WnormDegF);
    
    WnormDegS = mod(hdgS-outTableSlow.SwellWaveDir,360);
    hdgMinusWDirS=min(360-WnormDegS, WnormDegS);
    
    close all
    f2 = figure;
    set(gcf,'Position',[200 500 800 600]);
    hold on
    plot(hdgMinusWDirF,outTableFast.Std, 'o', 'MarkerFaceColor', 'b');
    plot(hdgMinusWDirS,outTableSlow.Std, 'oc', 'MarkerFaceColor', 'c');
    set(gca,'XLim',[0 180]);
    set(gca,'YLim',[0,1.5]);
    xlabel('abs(heading - swell wave dir) [deg]');
    ylabel('Liebe standard deviation [dB]');
    title([projectNameA ' Liebe standard deviation vs heading minus swell wave direction']);
    legend('Fast','Slow');
    plabel = [directories.figdir projectNameA '_biasStdVSwsbrdir_all'];
    set(gcf,'PaperPositionMode','auto')
    print(plabel,'-dpng','-r0')
    
end