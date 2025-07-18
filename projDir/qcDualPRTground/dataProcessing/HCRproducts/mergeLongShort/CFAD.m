% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

outdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

figdir=[indir(1:end-14),'mergeLongShort/cfad/'];

%% Run processing

% Go through iops
for ii=1:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];
    data.DBZ_short=[];
    data.FLAG_short=[];
    data.DBZ_long=[];
    data.FLAG_long=[];
           
    %% Load data
    disp('Loading data ...');

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.DBZ_short(data.FLAG_short~=1)=nan;
    data.DBZ_long(data.FLAG_long~=1)=nan;

    %% Count dbz values
    disp('Counting DBZ values ...');
    dbzEdges=-100:0.5:100;

    countsIOPs=zeros(size(data.range,1),length(dbzEdges)-1);
    countsIOPl=zeros(size(data.range,1),length(dbzEdges)-1);

    for jj=1:size(data.range,1)
        countsIOPs(jj,:)=histcounts(data.DBZ_short(jj,:),dbzEdges);
        countsIOPl(jj,:)=histcounts(data.DBZ_long(jj,:),dbzEdges);
    end

    if ii==1
        countsIOPallS=countsIOPs;
        countsIOPallL=countsIOPl;
    else
        countsIOPallS=countsIOPallS+countsIOPs;
        countsIOPallL=countsIOPallL+countsIOPl;
    end

end

countsIOPallS(countsIOPallS==0)=nan;
countsIOPallL(countsIOPallL==0)=nan;

%% Plot all iops together

% Short

close all
f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12);
t = tiledlayout(2,5,'TileSpacing','tight','Padding','compact');
colormap('jet');

surf(dbzEdges(1:end-1),data.range(:,1)./1000,countsIOPallS,'edgecolor','none');
view(2);

clim([0,20000])

xlim([-50,20])

colorbar

xlabel('Reflectivity (dBZ)')
ylabel('Range (km)')

box on
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_cfad_short'],'-dpng','-r0')

% Long
f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12);
t = tiledlayout(2,5,'TileSpacing','tight','Padding','compact');
colormap('jet');

surf(dbzEdges(1:end-1),data.range(:,1)./1000,countsIOPallL,'edgecolor','none');
view(2);

clim([0,20000])

xlim([-50,20])

colorbar

xlabel('Reflectivity (dBZ)')
ylabel('Range (km)')

box on
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_cfad_long'],'-dpng','-r0')

range=data.range;
save([outdir,'cfad.mat'],'dbzEdges','range','countsIOPallL');