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

%% Process regions
minArea=0.25;

%% Area of updraft regions

close all

edges=[0:0.25:5,inf];
xlab='Updraft area (km^2)';
figname=[figdir,'upRegs/upRegArea_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'area',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegArea_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'area',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);

%% Width of updraft regions

close all

edges=[0:1:14,inf];
xlab='Width of updraft regions (km)';
figname=[figdir,'upRegs/upRegWidth_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'width',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegWidth_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'width',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);

%% Depth of updraft regions

close all

edges=[0:0.2:2.8,inf];
xlab='Depth of updraft regions (km)';
figname=[figdir,'upRegs/upRegDepth_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'depth',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegDepth_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'depth',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean vel of updraft regions

close all

edges=[0:0.2:2.8,inf];
xlab='Mean velocity of updraft regions (m s^{-1})';
figname=[figdir,'upRegs/upRegMeanVel_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'meanVel',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegMeanVel_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'meanVel',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);

%% Max vel of updraft regions

close all

edges=[0:0.5:4.5,inf];
xlab='Max velocity of updraft regions (m s^{-1})';
figname=[figdir,'upRegs/upRegMaxVel_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'maxVel',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegMaxVel_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'maxVel',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);

%% Altitude of updraft regions

close all

edges=[0:10:100];
xlab='Altitude percentile of updraft regions (%)';
figname=[figdir,'upRegs/upRegAslPerc_',num2str(minArea),'.png'];
plotStatsRegsProjects(plotV.upRegsAll,'cloudAltPerc',minArea,edges,xlab,figname,classTypes,colmapCC);

figname=[figdir,'upRegs/upRegAslPerc_'];
plotStatsLocs_upRegs(plotV.upRegsAll,'cloudAltPerc',minArea,xlab,lonLims,latLims,figname,classTypes,colmapCC);
