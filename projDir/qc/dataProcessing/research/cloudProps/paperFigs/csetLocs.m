% Plot cset locations

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'cset'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));


qcVersion='v3.0';
quality='qc3';
ii=1;

cfDir=HCRdir(project{ii},quality,qcVersion,freqData);

indir=[cfDir(1:end-5),'cloudProps/'];

in.(project{ii})=load([indir,project{ii},'_cloudProps_2.mat']);

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
    'ConvYoungShallow','ConvYoungMid','ConvYoungDeep',...
    'ConvMatureShallow','ConvMatureMid','ConvMatureDeep'};

classTypeNames={'(a) CloudLow','(b) CloudMid','(e) CloudHigh',...
    '(b) StratShallow','(k) StratMid','(l) StratDeep',...
    '(d) ConvShallow','(e) ConvMid','(f) ConvDeep',...
    '(c) ConvStratShallow','(h) ConvStratMid','(i) ConvStratDeep'};

close all

%% Plot locations

close all

lonLims=[-160,-120;
    130,165;
    -95,-75];

latLims=[15,45;
    -65,-40
    -0,15];


lons=plotV.lonWholeAll;
lats=plotV.latWholeAll;

load coastlines

%% Per flight hour
load('/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/flightHourGrids.mat');

fig=figure('DefaultAxesFontSize',11,'position',[100,10,600,1270],'renderer','painters','visible','on');
colmap=turbo(105);
colormap([[0.8,0.8,0.8];colmap(1:100,:)]);

% Loop through projects
projects=fields(lons);
gridStep=1;

bb=1;

lonLength=(lonLims(bb,2)-lonLims(bb,1))/gridStep;
latLength=(latLims(bb,2)-latLims(bb,1))/gridStep;

lonSteps=lonLims(bb,1):gridStep:lonLims(bb,2);
latSteps=latLims(bb,1):gridStep:latLims(bb,2);

plotNum=num2str(1:12);

nums=[1,4,10,7,3];

for aa=1:length(nums)
    hourGrid=zeros(latLength,lonLength);
    s.plotNum{aa}=subplot(5,1,aa);

    hold on

    thisLons=lons.cset.(classTypes{nums(aa)});
    thisLats=lats.cset.(classTypes{nums(aa)});

    % Loop through output grid
    for ii=1:size(hourGrid,1)
        for jj=1:size(hourGrid,2)
            pixInds=find(thisLats>latSteps(ii) & thisLats<=latSteps(ii+1) & ...
                thisLons>lonSteps(jj) & thisLons<=lonSteps(jj+1));
            hzPix=length(pixInds);
            hourGrid(ii,jj)=hourGrid(ii,jj)+hzPix;
        end
    end

    perHour=hourGrid./flightHourGrids.cset;
    %perHour(perHour==0)=nan;
    perHour(flightHourGrids.cset<0.167)=nan;
    perHour=perHour./sum(perHour(:),'omitnan');

    h=imagesc(lonSteps(1:end-1)+gridStep/2,latSteps(1:end-1)+gridStep/2,perHour);
    set(h,'AlphaData',~isnan(perHour));

    caxis([-0.0009,0.1]);
    
    if aa==5
        cb=colorbar;
    end

    ylabel('Latitude (deg)');
    yticks(15:5:45);
    %yticklabels({'-95','-90','-85','-80'});

    if aa<5
        s.plotNum{aa}.XTickLabel=[];
    else
        xlabel('Longitude (deg)');
%         xticks(-95:5:-80);
%         xticklabels({'-95','-90','-85','-80'});
    end

    xlim([-160,-120]);
    ylim([18,44]);

    plot(coastlon,coastlat,'-k')
    text(-159,42,[classTypeNames{nums(aa)}],'fontsize',12,'FontWeight','bold','BackgroundColor','w');

    grid on
    box on
end

s.plotNum{1}.Position=[0.1 0.806 0.8 0.185];
s.plotNum{2}.Position=[0.1 0.614 0.8 0.185];
s.plotNum{3}.Position=[0.1 0.422 0.8 0.185];
s.plotNum{4}.Position=[0.1 0.23 0.8 0.185];
s.plotNum{5}.Position=[0.1 0.038 0.8 0.185];

cb.Position=[0.91,0.037,0.027,0.185];

set(gcf,'PaperPositionMode','auto')
print([figdir,'cset_locs_2.png'],'-dpng','-r0');
