% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

projects={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

for ii=1:length(projects)
    if strcmp(projects{ii},'cset')
        qcVersion='v3.0';
        quality='qc3';
    else
        qcVersion='v3.1';
        quality='qc3';
    end

    cfDir=HCRdir(projects{ii},quality,qcVersion,freqData);

    indir=[cfDir(1:end-5),'cloudProps/'];

    in.(projects{ii})=load([indir,projects{ii},'_cloudProps.mat']);
end

plotVars=fields(in.cset);

for ii=1:length(plotVars)
    for jj=1:length(projects)
        if ~(strcmp(projects{jj},'spicule') & strcmp(plotVars{ii},'sstAll'))
            plotV.(plotVars{ii}).(projects{jj})=in.(projects{jj}).(plotVars{ii});
        end
    end
end

%% Plot

disp('Plotting ...')


%% Plot

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,1100],'visible','on','renderer','painters');

load coastlines

% (a) Cloud base

s1=subplot(3,1,1);

hold on

thisLons=plotV.lonAll.cset.CloudLow;

thisLats=plotV.latAll.cset.CloudLow;

plotVar=plotV.cloudBaseAll.cset.CloudLow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 3]);
s1.Colormap=jet;
cb1=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(a) Cloud base (km). Low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (b) Cloud top

s2=subplot(3,1,2);

hold on

thisLons=cat(1,plotV.lonAll.cset.CloudLow,plotV.lonAll.cset.ConvYoungShallow, ...
    plotV.lonAll.cset.ConvMatureShallow,plotV.lonAll.cset.StratShallow);

thisLats=cat(1,plotV.latAll.cset.CloudLow,plotV.latAll.cset.ConvYoungShallow, ...
    plotV.latAll.cset.ConvMatureShallow,plotV.latAll.cset.StratShallow);

plotVar=cat(1,plotV.cloudTopAll.cset.CloudLow,plotV.cloudTopAll.cset.ConvYoungShallow, ...
    plotV.cloudTopAll.cset.ConvMatureShallow,plotV.cloudTopAll.cset.StratShallow);

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 4]);
s2.Colormap=jet;
cb2=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(b) Cloud top (km). Shallow and low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (c) Cloud width

s3=subplot(3,1,3);

hold on

thisLons=cat(1,plotV.lonAll.cset.CloudLow,plotV.lonAll.cset.ConvYoungShallow, ...
    plotV.lonAll.cset.ConvMatureShallow,plotV.lonAll.cset.StratShallow);

thisLats=cat(1,plotV.latAll.cset.CloudLow,plotV.latAll.cset.ConvYoungShallow, ...
    plotV.latAll.cset.ConvMatureShallow,plotV.latAll.cset.StratShallow);

plotVar=cat(1,plotV.cloudLengthAll.cset.CloudLow,plotV.cloudLengthAll.cset.ConvYoungShallow, ...
    plotV.cloudLengthAll.cset.ConvMatureShallow,plotV.cloudLengthAll.cset.StratShallow);

% thisLons(plotVar<20)=[];
% thisLats(plotVar<20)=[];
% plotVar(plotVar<20)=[];

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 80]);
s3.Colormap=jet;
cb3=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(c) Cloud width (km). Shallow and low clouds.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.095,0.71,0.8,0.27];
s2.Position=[0.095,0.377,0.8,0.27];
s3.Position=[0.095,0.045,0.8,0.27];

cb1.Position=[0.915,0.71,0.035,0.27];
cb2.Position=[0.915,0.377,0.035,0.27];
cb3.Position=[0.915,0.045,0.035,0.27];

set(gcf,'PaperPositionMode','auto')
print([figdir,'csetDims.png'],'-dpng','-r0');
