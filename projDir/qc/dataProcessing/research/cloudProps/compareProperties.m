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

%% Max refl

close all

edges=-50:5:30;
xlab='Maximum reflectivity (dBZ)';
figname=[figdir,'propsAll/maxRefl.png'];

plotStatsProjects(plotV.maxReflAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean refl

close all

edges=-50:5:30;
xlab='Mean reflectivity (dBZ)';
figname=[figdir,'propsAll/meanRefl.png'];

plotStatsProjects(plotV.meanReflAll,edges,xlab,figname,classTypes,colmapCC);

%% Max convectivity

close all

edges=0:0.1:1;
xlab='Max convectivity';
figname=[figdir,'propsAll/maxConv.png'];

plotStatsProjects(plotV.maxConvAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean convectivity

close all

edges=0:0.1:1;
xlab='Mean convectivity';
figname=[figdir,'propsAll/meanConv.png'];

plotStatsProjects(plotV.meanConvAll,edges,xlab,figname,classTypes,colmapCC);

%% Cloud depth

close all

edges=0:1:15;
xlab='Cloud depth (km)';
figname=[figdir,'propsAll/cloudDepth.png'];

plotStatsProjects(plotV.cloudDepthAll,edges,xlab,figname,classTypes,colmapCC);

%% Max temperature

close all

edges=-80:10:30;
xlab='Max temperature (C)';
figname=[figdir,'environment/maxTemp.png'];

plotStatsProjects(plotV.maxTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Min temperature

close all

edges=-80:10:30;
xlab='Min temperature (C)';
figname=[figdir,'environment/minTemp.png'];

plotStatsProjects(plotV.minTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean temperature

close all

edges=-80:10:30;
xlab='Mean temperature (C)';
figname=[figdir,'environment/meanTemp.png'];

plotStatsProjects(plotV.meanTempAll,edges,xlab,figname,classTypes,colmapCC);

%% Max pressure

close all

edges=100:100:1100;
xlab='Max pressure (hPa)';
figname=[figdir,'environment/maxPress.png'];

plotStatsProjects(plotV.maxPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Min pressure

close all

edges=100:100:1100;
xlab='Min pressure (hPa)';
figname=[figdir,'environment/minPress.png'];

plotStatsProjects(plotV.minPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean pressure

close all

edges=100:100:1100;
xlab='Mean pressure (hPa)';
figname=[figdir,'environment/meanPress.png'];

plotStatsProjects(plotV.meanPressAll,edges,xlab,figname,classTypes,colmapCC);

%% Icing level

close all

edges=0:0.5:6;
xlab='Icing level (km)';
figname=[figdir,'environment/iceLev.png'];

plotStatsProjects(plotV.iceLevAll,edges,xlab,figname,classTypes,colmapCC);

%% SST

close all

edges=-5:2:35;
xlab='SST (C)';
figname=[figdir,'environment/sst.png'];

plotStatsProjects(plotV.sstAll,edges,xlab,figname,classTypes,colmapCC);

%% Updraft fraction

close all

edges=0:0.1:1;
xlab='Updraft fraction';
figname=[figdir,'propsAll/upFrac.png'];

plotStatsProjects(plotV.upFracAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum strength of updraft

close all

edges=[0:1:14,inf];
xlab='Max up velocity (m s^{-1})';
figname=[figdir,'propsAll/upMaxStrength.png'];

plotStatsProjects(plotV.upMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);

%% Maximum strength of downdraft

close all

edges=[0:1:14,inf];
xlab='Max down velocity (m s^{-1})';
figname=[figdir,'propsAll/downMaxStrength.png'];

plotStatsProjects(plotV.downMaxStrengthAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean strength of updraft

close all

edges=[0:0.2:3.8,inf];
xlab='Mean up velocity (m s^{-1})';
figname=[figdir,'propsAll/upMeanStrength.png'];

plotStatsProjects(plotV.upMeanStrengthAll,edges,xlab,figname,classTypes,colmapCC);

%% Mean strength of downdraft

close all

edges=[0:0.2:3.8,inf];
xlab='Mean down velocity (m s^{-1})';
figname=[figdir,'propsAll/downMeanStrength.png'];

plotStatsProjects(plotV.downMeanStrengthAll,edges,xlab,figname,classTypes,colmapCC);

%% Process regions
minPixNum=1000;

%% Pixel number updraft regions

close all

edges=[0:1000:9000,inf];
xlab='Pixels in updraft regions';
figname=[figdir,'upRegs/upRegPixNum_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'numPix',minPixNum,edges,xlab,figname,classTypes,colmapCC);

%% Width of updraft regions

close all

edges=[0:1:14,inf];
xlab='Width of updraft regions (km)';
figname=[figdir,'upRegs/upRegWidth_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'width',minPixNum,edges,xlab,figname,classTypes,colmapCC);

%% Depth of updraft regions

close all

edges=[0:0.2:2.8,inf];
xlab='Depth of updraft regions (km)';
figname=[figdir,'upRegs/upRegDepth_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'depth',minPixNum,edges,xlab,figname,classTypes,colmapCC);

%% Mean vel of updraft regions

close all

edges=[0:0.2:2.8,inf];
xlab='Mean velocity of updraft regions (m s^{-1})';
figname=[figdir,'upRegs/upRegMeanVel_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'meanVel',minPixNum,edges,xlab,figname,classTypes,colmapCC);

%% Max vel of updraft regions

close all

edges=[0:0.5:4.5,inf];
xlab='Max velocity of updraft regions (m s^{-1})';
figname=[figdir,'upRegs/upRegMaxVel_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'maxVel',minPixNum,edges,xlab,figname,classTypes,colmapCC);

%% Altitude of updraft regions

close all

edges=[0:10:100];
xlab='Altitude percentile of updraft regions (%)';
figname=[figdir,'upRegs/upRegAslPerc_',num2str(minPixNum),'.png'];

plotStatsRegsProjects(plotV.upRegsAll,'cloudAltPerc',minPixNum,edges,xlab,figname,classTypes,colmapCC);

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