% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'cset','socrates','otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

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

    in.(project{ii})=load([indir,project{ii},'_ssts.mat']);
end


%% Plot

disp('Plotting ...')

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];

load coastlines

fig=figure('DefaultAxesFontSize',11,'position',[100,100,600,1100],'renderer','painters','visible','on');
colormap('jet');

% CSET
s1=subplot(3,1,1);
hold on

inds=1:100:length(in.cset.latAll);
scatter(in.cset.lonAll(inds),in.cset.latAll(inds),1,in.cset.sstAll(inds),'filled');

inds2=find(isnan(in.cset.sstAll));
scatter(in.cset.lonAll(inds2),in.cset.latAll(inds2),1,[0.5,0.5,0.5],'filled');

caxis([15 30]);

xlim(lonLims(1,:));
ylim(latLims(1,:));

plot(coastlon,coastlat,'-k')

grid on
box on

cb1=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(a) CSET, SST (',char(176),'C)']);

% SOCRATES

s2=subplot(3,1,2);
hold on

inds=1:100:length(in.socrates.latAll);
scatter(in.socrates.lonAll(inds),in.socrates.latAll(inds),1,in.socrates.sstAll(inds),'filled');

inds2=find(isnan(in.socrates.sstAll));
scatter(in.socrates.lonAll(inds2),in.socrates.latAll(inds2),1,[0.5,0.5,0.5],'filled');

caxis([0 20]);

xlim(lonLims(2,:));
ylim(latLims(2,:));

plot(coastlon,coastlat,'-k')

grid on
box on

cb2=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(b) SOCRATES, SST (',char(176),'C)']);

% OTREC

s3=subplot(3,1,3);
hold on

inds=1:100:length(in.otrec.latAll);
scatter(in.otrec.lonAll(inds),in.otrec.latAll(inds),1,in.otrec.sstAll(inds),'filled');

inds2=find(isnan(in.otrec.sstAll));
scatter(in.otrec.lonAll(inds2),in.otrec.latAll(inds2),1,[0.5,0.5,0.5],'filled');

caxis([25 30]);

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')

grid on
box on

cb3=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(c) OTREC, SST (',char(176),'C)']);

s1.Position=[0.095,0.71,0.8,0.27];
s2.Position=[0.095,0.377,0.8,0.27];
s3.Position=[0.095,0.045,0.8,0.27];

cb1.Position=[0.915,0.71,0.035,0.27];
cb2.Position=[0.915,0.377,0.035,0.27];
cb3.Position=[0.915,0.045,0.035,0.27];

set(gcf,'PaperPositionMode','auto')
print([figdir,'sstTracks.png'],'-dpng','-r0');
