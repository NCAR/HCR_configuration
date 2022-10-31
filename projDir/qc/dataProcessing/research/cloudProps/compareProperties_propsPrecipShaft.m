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

%% Shafts total length

minLength=2.5;

close all

edges=[2.5:10:100];
xlab='Total length of precip shafts (km)';
figname=[figdir,'precShafts/precShafts_totLength'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'shaftKM',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'shaftKM',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_totLength_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'shaftKM',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Precip fraction

close all

edges=[0:0.1:1];
xlab='Precip fraction';
figname=[figdir,'precShafts/precShafts_frac'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'frac',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'frac',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_frac_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'frac',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean refl

close all

edges=[-inf,-30:5:10,inf];
xlab='Mean precip reflectivity (dBZ)';
figname=[figdir,'precShafts/precShafts_meanRefl'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'meanRef',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'meanRef',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_meanRefl_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'meanRef',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Max refl

close all

edges=[-inf,-30:5:15,inf];
xlab='Max precip reflectivity (dBZ)';
figname=[figdir,'precShafts/precShafts_maxRefl'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'maxRefl',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'maxRefl',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_maxRefl_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'maxRefl',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Mean vel

close all

edges=[-inf,0:0.5:5,inf];
xlab='Mean precip velocity (m s^{-1})';
figname=[figdir,'precShafts/precShafts_meanVel'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'meanVel',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'meanVel',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_meanVel_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'meanVel',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);

%% Max vel

close all

edges=[-inf,0:1:10,inf];
xlab='Max precip velocity (m s^{-1})';
figname=[figdir,'precShafts/precShafts_maxVel'];
plotStatsPrecShaftsProjects(plotV.precShaftsAll,'maxVel',minLength,edges,xlab,figname,classTypes,colmapCC);
plotStatsPrecShaftsProjects_box(plotV.precShaftsAll,'maxVel',minLength,xlab,figname,classTypes,colmapCC,[edges(1),edges(end)]);

figname=[figdir,'precShafts/precShafts_maxVel_'];
plotStatsLocs_precShafts(plotV.precShaftsAll,'maxVel',minLength,xlab,plotV.lonAll,plotV.latAll,lonLims,latLims,figname,classTypes,colmapCC);
