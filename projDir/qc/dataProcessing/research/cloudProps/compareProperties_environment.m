% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

for ii=1:length(project)
    if strcmp(project{ii},'spicule')
        qcVersion='v1.1';
        quality='qc1';
    elseif strcmp(project{ii},'cset')
        qcVersion='v3.0';
        quality='qc3';
    elseif strcmp(project{ii},'noreaster')
        qcVersion='v2.0';
        quality='qc2';
    else
        qcVersion='v3.1';
        quality='qc3';
    end

    cfDir=HCRdir(project{ii},quality,qcVersion,freqData);

    indir=[cfDir(1:end-5),'cloudProps/'];

    in.(project{ii})=load([indir,project{ii},'_cloudProps.mat']);
end

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(project)
        if ~(strcmp(project{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(project{jj})=in.(project{jj}).(plotVars{ii});
        end
    end
end

%% Plot

disp('Plotting ...')

classTypes={'CloudLow','CloudMid','CloudHigh',...
        'StratShallow','StratMid','StratDeep',...
        'ConvYoungShallow','ConvYoungMid','ConvYongDeep',...
        'ConvMatureShallow','ConvMatureMid','ConvMatureDeep'};

colmapCC=[204,255,204;
    153,204,0;
    0,128,0;
    0,204,255;
    51,102,255;
    0,0,180;
    255,204,0;
    255,102,0;
    220,0,0;
    255,153,220;
    204,153,255;
    128,0,128];

colmapCC=colmapCC./255;

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

%% Max temperature

close all

edges=-80:10:30;
xlab='Max temperature (C)';
figname=[figdir,'environment/maxTemp.png'];
plotStatsProjects(plotV.maxTempAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/maxTemp_'];
plotStatsLocs(plotV.maxTempAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Min temperature

close all

edges=-80:10:30;
xlab='Min temperature (C)';
figname=[figdir,'environment/minTemp.png'];
plotStatsProjects(plotV.minTempAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/minTemp_'];
plotStatsLocs(plotV.minTempAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean temperature

close all

edges=-80:10:30;
xlab='Mean temperature (C)';
figname=[figdir,'environment/meanTemp.png'];

plotStatsProjects(plotV.meanTempAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/meanTemp_'];
plotStatsLocs(plotV.meanTempAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Max pressure

close all

edges=100:100:1100;
xlab='Max pressure (hPa)';
figname=[figdir,'environment/maxPress.png'];
plotStatsProjects(plotV.maxPressAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/maxPress_'];
plotStatsLocs(plotV.maxPressAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Min pressure

close all

edges=100:100:1100;
xlab='Min pressure (hPa)';
figname=[figdir,'environment/minPress.png'];
plotStatsProjects(plotV.minPressAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/minPress_'];
plotStatsLocs(plotV.minPressAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean pressure

close all

edges=100:100:1100;
xlab='Mean pressure (hPa)';
figname=[figdir,'environment/meanPress.png'];
plotStatsProjects(plotV.meanPressAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/meanPress_'];
plotStatsLocs(plotV.meanPressAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Icing level

close all

edges=0:0.5:6;
xlab='Icing level (km)';
figname=[figdir,'environment/iceLev.png'];
plotStatsProjects(plotV.iceLevAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/iceLev_'];
plotStatsLocs(plotV.iceLevAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% SST

close all

edges=-5:2:35;
xlab='SST (C)';
figname=[figdir,'environment/sst.png'];
plotStatsProjects(plotV.sstAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'environment/sst_'];
plotStatsLocs(plotV.sstAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);
