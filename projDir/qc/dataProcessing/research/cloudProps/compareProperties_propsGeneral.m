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

close all

%% Plot locations

close all

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

figname=[figdir,'locations.png'];

plotLocsProjects(plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Plot types locations

close all

figname=[figdir,'locationsTypes_'];

plotLocsTypes(plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Max refl

close all

edges=-50:5:30;
xlab='Maximum reflectivity (dBZ)';
figname=[figdir,'propsAll/maxRefl.png'];
plotStatsProjects(plotV.maxReflAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/maxRefl_'];
plotStatsLocs(plotV.maxReflAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean refl

close all

edges=-50:5:30;
xlab='Mean reflectivity (dBZ)';
figname=[figdir,'propsAll/meanRefl.png'];
plotStatsProjects(plotV.meanReflAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/meanRefl_'];
plotStatsLocs(plotV.meanReflAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Max convectivity

close all

edges=0:0.1:1;
xlab='Max convectivity';
figname=[figdir,'propsAll/maxConv.png'];
plotStatsProjects(plotV.maxConvAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/maxConv_'];
plotStatsLocs(plotV.maxConvAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean convectivity

close all

edges=0:0.1:1;
xlab='Mean convectivity';
figname=[figdir,'propsAll/meanConv.png'];
plotStatsProjects(plotV.meanConvAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/meanConv_'];
plotStatsLocs(plotV.meanConvAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Cloud depth

close all

edges=0:1:15;
xlab='Cloud depth (km)';
figname=[figdir,'propsAll/cloudDepth.png'];
plotStatsProjects(plotV.cloudDepthAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/cloudDepth_'];
plotStatsLocs(plotV.cloudDepthAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Cloud top

close all

edges=0:1:15;
xlab='Cloud top (km)';
figname=[figdir,'propsAll/cloudTop.png'];
plotStatsProjects(plotV.cloudTopAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/cloudTop_'];
plotStatsLocs(plotV.cloudTopAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Cloud base

close all

edges=0:1:15;
xlab='Cloud base (km)';
figname=[figdir,'propsAll/cloudBase.png'];
plotStatsProjects(plotV.cloudBaseAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/cloudBase_'];
plotStatsLocs(plotV.cloudBaseAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Cloud layers

close all

edges=0:1:10;
xlab='Cloud layers';
figname=[figdir,'propsAll/cloudLayers.png'];
plotStatsProjects(plotV.cloudLayersAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/cloudLayers_'];
plotStatsLocs(plotV.cloudLayersAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Updraft fraction

close all

edges=0:0.1:1;
xlab='Updraft fraction';
figname=[figdir,'propsAll/upFrac.png'];
plotStatsProjects(plotV.upFracAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/upFrac_'];
plotStatsLocs(plotV.upFracAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Maximum strength of updraft

close all

edges=[0:1:14,inf];
xlab='Max up velocity (m s^{-1})';
figname=[figdir,'propsAll/upMaxStrength.png'];
plotStatsProjects(plotV.upMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/upMaxStrength_'];
plotStatsLocs(plotV.upMaxStrengthAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Maximum strength of downdraft

close all

edges=[0:1:14,inf];
xlab='Max down velocity (m s^{-1})';
figname=[figdir,'propsAll/downMaxStrength.png'];
plotStatsProjects(plotV.downMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/downMaxStrength_'];
plotStatsLocs(plotV.downMaxStrengthAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean strength of updraft

close all

edges=[0:0.2:3.8,inf];
xlab='Mean up velocity (m s^{-1})';
figname=[figdir,'propsAll/upMeanStrength.png'];
plotStatsProjects(plotV.upMeanStrengthAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/upMeanStrength_'];
plotStatsLocs(plotV.upMeanStrengthAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean strength of downdraft

close all

edges=[0:0.2:3.8,inf];
xlab='Mean down velocity (m s^{-1})';
figname=[figdir,'propsAll/downMeanStrength.png'];
plotStatsProjects(plotV.downMeanStrengthAll,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'propsAll/downMeanStrength_'];
plotStatsLocs(plotV.downMeanStrengthAll,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);
