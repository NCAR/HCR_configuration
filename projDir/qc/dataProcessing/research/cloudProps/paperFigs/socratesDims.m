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

plotVars=fields(in.socrates);

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

colm=jet(64);
colm=colm(1:55,:);

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

close all


fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,700],'visible','on','renderer','painters');

load coastlines

% (a) Cloud top

s1=subplot(2,1,1);

hold on

thisLons=cat(1,plotV.lonAll.socrates.CloudLow,plotV.lonAll.socrates.ConvYoungShallow, ...
    plotV.lonAll.socrates.ConvMatureShallow,plotV.lonAll.socrates.StratShallow,...
    plotV.lonAll.socrates.ConvYoungMid,plotV.lonAll.socrates.StratMid);

thisLats=cat(1,plotV.latAll.socrates.CloudLow,plotV.latAll.socrates.ConvYoungShallow, ...
    plotV.latAll.socrates.ConvMatureShallow,plotV.latAll.socrates.StratShallow,...
    plotV.latAll.socrates.ConvYoungMid,plotV.latAll.socrates.StratMid);

plotVar=cat(1,plotV.cloudTopAll.socrates.CloudLow,plotV.cloudTopAll.socrates.ConvYoungShallow, ...
    plotV.cloudTopAll.socrates.ConvMatureShallow,plotV.cloudTopAll.socrates.StratShallow,...
    plotV.cloudTopAll.socrates.ConvYoungMid,plotV.cloudTopAll.socrates.StratMid);

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 3.5]);
s1.Colormap=colm;
cb1=colorbar;

xlim(lonLims(2,:));
ylim(latLims(2,:));

plot(coastlon,coastlat,'-k')
title('(a) Cloud top (km). Precipitating Shallow, Mid, and CloudLow.');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

% (c) Cloud width

s2=subplot(2,1,2);

hold on

thisLons=cat(1,plotV.lonAll.socrates.CloudLow,plotV.lonAll.socrates.ConvYoungShallow, ...
    plotV.lonAll.socrates.CloudMid,plotV.lonAll.socrates.CloudHigh,...
    plotV.lonAll.socrates.ConvMatureShallow,plotV.lonAll.socrates.StratShallow,...
    plotV.lonAll.socrates.ConvYoungMid,plotV.lonAll.socrates.StratMid);

thisLats=cat(1,plotV.latAll.socrates.CloudLow,plotV.latAll.socrates.ConvYoungShallow, ...
    plotV.latAll.socrates.CloudMid,plotV.latAll.socrates.CloudHigh,...
    plotV.latAll.socrates.ConvMatureShallow,plotV.latAll.socrates.StratShallow,...
    plotV.latAll.socrates.ConvYoungMid,plotV.latAll.socrates.StratMid);

plotVar=cat(1,plotV.cloudLengthAll.socrates.CloudLow,plotV.cloudLengthAll.socrates.ConvYoungShallow, ...
    plotV.cloudLengthAll.socrates.CloudMid,plotV.cloudLengthAll.socrates.CloudHigh,...
    plotV.cloudLengthAll.socrates.ConvMatureShallow,plotV.cloudLengthAll.socrates.StratShallow,...
    plotV.cloudLengthAll.socrates.ConvYoungMid,plotV.cloudLengthAll.socrates.StratMid);

% thisLons(plotVar<20)=[];
% thisLats(plotVar<20)=[];
% plotVar(plotVar<20)=[];

scatter(thisLons,thisLats,25,plotVar,'filled');

caxis([0 80]);
s2.Colormap=colm;
cb2=colorbar;

xlim(lonLims(2,:));
ylim(latLims(2,:));

plot(coastlon,coastlat,'-k')
title('(b) Cloud width (km).');

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

grid on
box on

s1.Position=[0.095,0.565,0.8,0.4];
s2.Position=[0.095,0.07,0.8,0.4];

cb1.Position=[0.915,0.565,0.035,0.4];
cb2.Position=[0.915,0.07,0.035,0.4];

set(gcf,'PaperPositionMode','auto')
print([figdir,'socratesDims.png'],'-dpng','-r0');
