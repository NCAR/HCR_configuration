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

close all

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

load coastlines

colmap=jet;
colmapC=cat(1,[0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0; ...
    0,0,0.05;0,0,0.1;0,0,0.15;0,0,0.2;0,0,0.25;0,0,0.3;0,0,0.35;0,0,0.4;0,0,0.45;0,0,0.5],colmap(1:230,:));

fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,700],'visible','on','renderer','painters');

% (a) CSET

s1=subplot(2,1,1);

hold on

thisLons=[plotV.lonAll.cset.StratShallow;plotV.lonAll.cset.ConvYoungShallow;plotV.lonAll.cset.ConvMatureShallow];
thisLats=[plotV.latAll.cset.StratShallow;plotV.latAll.cset.ConvYoungShallow;plotV.latAll.cset.ConvMatureShallow];

plotVar=[plotV.precShaftsAll.cset.StratShallow.maxVel;plotV.precShaftsAll.cset.ConvYoungShallow.maxVel;plotV.precShaftsAll.cset.ConvMatureShallow.maxVel];

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 7]);
s1.Colormap=colmapC;
cb1=colorbar;

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')
title('(a) CSET maximum velocity (m s^{-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (b) OTREC

s2=subplot(2,1,2);

hold on

thisLons=[plotV.lonAll.otrec.StratShallow;plotV.lonAll.otrec.ConvYoungShallow;plotV.lonAll.otrec.ConvMatureShallow];
thisLats=[plotV.latAll.otrec.StratShallow;plotV.latAll.otrec.ConvYoungShallow;plotV.latAll.otrec.ConvMatureShallow];

plotVar=[plotV.precShaftsAll.otrec.StratShallow.maxVel;plotV.precShaftsAll.otrec.ConvYoungShallow.maxVel;plotV.precShaftsAll.otrec.ConvMatureShallow.maxVel];

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 7.5]);
s2.Colormap=colmapC;
cb2=colorbar;

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')
title('(b) OTREC maximum velocity (m s{^-1})');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.095,0.565,0.8,0.4];
s2.Position=[0.095,0.07,0.8,0.4];

cb1.Position=[0.91,0.565,0.035,0.4];
cb2.Position=[0.91,0.07,0.035,0.4];

set(gcf,'PaperPositionMode','auto')
print([figdir,'shaftVel.png'],'-dpng','-r0');
