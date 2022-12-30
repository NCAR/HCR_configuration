% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

projects={'cset';'socrates';'otrec'}; %socrates, aristo, cset, otrec
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

cols3=[0,1,0; ...
    0,0,1; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,0,0; ...
    0,1,0; ...
    0,0,1; ...
    1,0,0];

projects3=cat(1,projects,projects,projects);

close all

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,800],'visible','on','renderer','painters');

% (a) Mean convectivity conv young shallow CSET

s1=subplot(2,2,1);

hold on

thisLons=plotV.lonAll.cset.ConvYoungShallow;
thisLats=plotV.latAll.cset.ConvYoungShallow;

plotVar=plotV.meanConvAll.cset.ConvYoungShallow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 1]);
colmap=jet;
colmapC=cat(1,[0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0.05;0,0,0.1;0,0,0.15;0,0,0.2;0,0,0.25;0,0,0.3;0,0,0.35;0,0,0.4;0,0,0.45;0,0,0.5],colmap(1:230,:));
s1.Colormap=colmapC;
cb1=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(a) CSET: mean convectivity');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (b) Mean convectivity conv young shallow OTREC

s2=subplot(2,2,2);

hold on

thisLons=plotV.lonAll.otrec.ConvYoungShallow;
thisLats=plotV.latAll.otrec.ConvYoungShallow;

plotVar=plotV.meanConvAll.otrec.ConvYoungShallow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 1]);
s2.Colormap=colmapC;
cb2=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(b) OTREC: mean convectivity');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (c) Mean up vel conv young shallow CSET

s3=subplot(2,2,3);

hold on

thisLons=plotV.lonAll.cset.ConvYoungShallow;
thisLats=plotV.latAll.cset.ConvYoungShallow;

plotVar=plotV.upMeanStrengthAll.cset.ConvYoungShallow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 1]);
colmap=jet;
s3.Colormap=colmap(10:230,:);
cb3=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(c) CSET: mean upward motion (m s^{-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (d) Mean up vel conv young shallow OTREC

s4=subplot(2,2,4);

hold on

thisLons=plotV.lonAll.otrec.ConvYoungShallow;
thisLats=plotV.latAll.otrec.ConvYoungShallow;

plotVar=plotV.upMeanStrengthAll.otrec.ConvYoungShallow;

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 1.2]);
colmap=jet;
s4.Colormap=colmap(10:230,:);
cb4=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(d) OTREC: mean upward motion (m s^{-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.045,0.56,0.4,0.4];
s2.Position=[0.545,0.56,0.4,0.4];
s3.Position=[0.045,0.065,0.4,0.4];
s4.Position=[0.545,0.065,0.4,0.4];

cb1.Position=[0.455,0.56,0.02,0.4];
cb2.Position=[0.955,0.56,0.02,0.4];
cb3.Position=[0.455,0.065,0.02,0.4];
cb4.Position=[0.955,0.065,0.02,0.4];

set(gcf,'PaperPositionMode','auto')
print([figdir,'convSpatial.png'],'-dpng','-r0');
