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

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,900],'renderer','painters','visible','on');
colormap('jet');

% WORLD
s1=subplot(2,2,1);
hold on

xlim([-180,180]);
ylim([-90,90]);

colsAll=[0,1,0;0,0,1;1,0,0];

plot(coastlon,coastlat,'-k')
for aa=1:3
    plot(lonLims(aa,:),[latLims(aa,1),latLims(aa,1)],'-','Color',colsAll(aa,:),'LineWidth',2);
    plot(lonLims(aa,:),[latLims(aa,2),latLims(aa,2)],'-','Color',colsAll(aa,:),'LineWidth',2);
    plot([lonLims(aa,1),lonLims(aa,1)],latLims(aa,:),'-','Color',colsAll(aa,:),'LineWidth',2);
    plot([lonLims(aa,2),lonLims(aa,2)],latLims(aa,:),'-','Color',colsAll(aa,:),'LineWidth',2);
end

text(lonLims(1,1),latLims(1,1)-8,'CSET','Color','g','FontSize',12,'FontWeight','bold');
text(lonLims(2,1)-40,latLims(2,1)-8,'SOCRATES','Color','b','FontSize',12,'FontWeight','bold');
text(lonLims(3,1)-20,latLims(3,1)-8,'OTREC','Color','r','FontSize',12,'FontWeight','bold');

grid on
box on

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(a) Project locations']);

% CSET
s2=subplot(2,2,2);
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

cb2=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(b) CSET, SST (',char(176),'C)']);

% SOCRATES

s3=subplot(2,2,3);
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

cb3=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(c) SOCRATES, SST (',char(176),'C)']);

% OTREC

s4=subplot(2,2,4);
hold on

inds=1:100:length(in.otrec.latAll);
scatter(in.otrec.lonAll(inds),in.otrec.latAll(inds),1,in.otrec.sstAll(inds),'filled');

inds2=find(isnan(in.otrec.sstAll));
scatter(in.otrec.lonAll(inds2),in.otrec.latAll(inds2),1,[0.5,0.5,0.5],'filled');

caxis([25 30]);

xlim(lonLims(3,:));
ylim(latLims(3,:));

plot(coastlon,coastlat,'-k')

plot([-94.8,-84.5],[2.5,2.5],'-k','LineWidth',2);
plot([-94.8,-84.5],[13.5,13.5],'-k','LineWidth',2);
plot([-94.8,-94.8],[2.5,13.5],'-k','LineWidth',2);
plot([-84.5,-84.5],[2.5,13.5],'-k','LineWidth',2);

plot([-81.7,-77.5],[2.5,2.5],'-k','LineWidth',2);
plot([-81.7,-77.5],[8,8],'-k','LineWidth',2);
plot([-81.7,-81.7],[2.5,8],'-k','LineWidth',2);
plot([-77.5,-77.5],[2.5,8],'-k','LineWidth',2);

plot([-83.4,-77.5],[9.1,9.1],'-k','LineWidth',2);
plot([-83.4,-77.5],[13.8,13.8],'-k','LineWidth',2);
plot([-77.5,-77.5],[9.1,13.8],'-k','LineWidth',2);
plot([-83.4,-83.4],[9.1,13.8],'-k','LineWidth',2);

text(-91.2,14.3,'ITCZ region','FontSize',12,'FontWeight','bold');
text(-83.9,14.5,'Caribbean region','FontSize',12,'FontWeight','bold');
text(-82.3,1.5,[{'Colombian'};{'coastal region'}],'FontSize',12,'FontWeight','bold');

grid on
box on

cb4=colorbar;

xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

title(['(d) OTREC, SST (',char(176),'C)']);

s1.Position=[0.05,0.56,0.43,0.41];
s2.Position=[0.54,0.56,0.38,0.41];
s3.Position=[0.05,0.06,0.38,0.41];
s4.Position=[0.54,0.06,0.38,0.41];

set(gcf,'PaperPositionMode','auto')
print([figdir,'sstTracks.png'],'-dpng','-r0');
