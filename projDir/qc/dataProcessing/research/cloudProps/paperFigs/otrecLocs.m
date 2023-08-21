% Plot OTREC locations

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project={'otrec'}; %socrates, aristo, cset, otrec
freqData='10hz';
whichModel='era5';

figdir='/scr/snow2/rsfdata/projects/cset/hcr/qc3/cfradial/v3.0_full/cloudPropsProjects/paperFigs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));


qcVersion='v3.1';
quality='qc3';
ii=1;

cfDir=HCRdir(project{ii},quality,qcVersion,freqData);

indir=[cfDir(1:end-5),'cloudProps/'];

in.(project{ii})=load([indir,project{ii},'_cloudProps_2.mat']);

plotVars=fields(in.otrec);

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

classTypeNames={'(a) CloudLow','(b) CloudMid','(c) CloudHigh',...
    '(j) StratShallow','(k) StratMid','(l) StratDeep',...
    '(d) ConvShallow','(e) ConvMid','(f) ConvDeep',...
    '(g) ConvStratShallow','(h) ConvStratMid','(i) ConvStratDeep'};

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

fig=figure('DefaultAxesFontSize',11,'position',[100,100,1200,1000],'renderer','painters','visible','on');
colmap=turbo(105);
colormap([[0.8,0.8,0.8];colmap(1:100,:)]);

% Loop through projects
projects=fields(lons);
gridStep=1;

bb=3;

lonLength=(lonLims(bb,2)-lonLims(bb,1))/gridStep;
latLength=(latLims(bb,2)-latLims(bb,1))/gridStep;

lonSteps=lonLims(bb,1):gridStep:lonLims(bb,2);
latSteps=latLims(bb,1):gridStep:latLims(bb,2);

plotNum=num2str(1:12);

nums=[1,2,3,7,8,9,10,11,12,4,5,6];

for aa=1:length(classTypes)
    hourGrid=zeros(latLength,lonLength);
    s.plotNum{aa}=subplot(4,3,aa);

    hold on

    thisLons=lons.otrec.(classTypes{nums(aa)});
    thisLats=lats.otrec.(classTypes{nums(aa)});

    % Loop through output grid
    for ii=1:size(hourGrid,1)
        for jj=1:size(hourGrid,2)
            pixInds=find(thisLats>latSteps(ii) & thisLats<=latSteps(ii+1) & ...
                thisLons>lonSteps(jj) & thisLons<=lonSteps(jj+1));
            hzPix=length(pixInds);
            hourGrid(ii,jj)=hourGrid(ii,jj)+hzPix;
        end
    end

    perHour=hourGrid./flightHourGrids.otrec;
    %perHour(perHour==0)=nan;
    perHour(flightHourGrids.otrec<0.5)=nan;
    perHour=perHour./sum(perHour(:),'omitnan');

    h=imagesc(lonSteps(1:end-1)+gridStep/2,latSteps(1:end-1)+gridStep/2,perHour);
    set(h,'AlphaData',~isnan(perHour));

    caxis([-0.0009,0.1]);
    if aa==12
        cb=colorbar;
    end

    if aa==2 | aa==3 | aa==5 | aa==6 | aa==8 | aa==9 | aa==11 | aa==12
        s.plotNum{aa}.YTickLabel=[];
    else
        ylabel('Latitude (deg)');
        %yticklabels({'-95','-90','-85','-80'});
    end
    if aa<10
        s.plotNum{aa}.XTickLabel=[];
    else
        xlabel('Longitude (deg)');
        xticks(-95:5:-80);
        xticklabels({'-95','-90','-85','-80'});
    end

    xlim(lonLims(bb,:));
    ylim(latLims(bb,:));

    plot(coastlon,coastlat,'-k')
    text(-94.5,13.5,[classTypeNames{nums(aa)}],'fontsize',12,'FontWeight','bold','BackgroundColor','w');

    grid on
    box on
end

s.plotNum{1}.Position=[0.042 0.765 0.295 0.225];
s.plotNum{2}.Position=[0.3475 0.765 0.295 0.225];
s.plotNum{3}.Position=[0.655 0.765 0.295 0.225];

s.plotNum{4}.Position=[0.042 0.525 0.295 0.225];
s.plotNum{5}.Position=[0.3475 0.525 0.295 0.225];
s.plotNum{6}.Position=[0.655 0.525 0.295 0.225];

s.plotNum{7}.Position=[0.042 0.285 0.295 0.225];
s.plotNum{8}.Position=[0.3475 0.285 0.295 0.225];
s.plotNum{9}.Position=[0.655 0.285 0.295 0.225];

s.plotNum{10}.Position=[0.042 0.045 0.295 0.225];
s.plotNum{11}.Position=[0.3475 0.045 0.295 0.225];
s.plotNum{12}.Position=[0.655 0.045 0.295 0.225];

cb.Position=[0.955,0.045,0.013,0.225];

set(gcf,'PaperPositionMode','auto')
print([figdir,'otrec_locs_2.png'],'-dpng','-r0');
