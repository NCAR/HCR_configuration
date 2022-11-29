% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

startTime=datetime(2019,9,21,16,30,0);
endTime=datetime(2019,9,21,16,42,0);

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=15;

saveFig=1;

figdir=[indir(1:end-5),'multScatt/cases/'];
if ~exist(figdir, 'dir')
    mkdir(figdir)
end

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/multScatt_',project,'.txt'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);

    %% Get data

    disp('Reading data ...');

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ=[];
    data.DBMVC=[];
    data.DBMHX=[];
    data.LDR=[];
    data.ICING_LEVEL=[];

    data=read_HCR(fileList,data,startTime,endTime);

    data.DBMHX(isnan(data.LDR))=-106;
    data.DBMHX(isnan(data.DBZ))=nan;

    data.DBMVC(isnan(data.DBZ))=nan;

    %% Find icing regions
    [LDRgrad,smoothLDR,iceMask]=findIce(data);

    %% Plot LDR, smooth and gradient
    close all
    disp('Plotting LDR ...');

    plotIndNum=1000;
    plotSpace=round(length(data.time)/plotIndNum);
    plotInds=1:plotSpace:length(data.time);

    fig=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',14);

    colormap('jet');

    s1=subplot(3,1,1);

    surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.LDR(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-25,0]);
    colorbar

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('LDR (dB)')

    s2=subplot(3,1,2);

    surf(data.time(plotInds),data.asl(:,plotInds)./1000,smoothLDR(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-25,0]);
    colorbar

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('Smoothed LDR (dB)')

    s3=subplot(3,1,3);
    hold on

    surf(data.time(plotInds),data.asl(1:end-1,plotInds)./1000,LDRgrad(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-0.2,0.4]);
    colorbar

    bounds=bwboundaries(iceMask);
    for kk=1:length(bounds)
        boundary=bounds{kk};
        btimes=data.time(boundary(:,2));
        linInds=sub2ind(size(data.LDR),boundary(:,1),boundary(:,2));
        basl=data.asl(linInds)./1000;
        plot(btimes,basl,'k','LineWidth',1.5);
    end

    s3.SortMethod='childorder';

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('LDR gradient (dB)')

    linkaxes([s1 s2 s3],'xy')

    if saveFig
        set(gcf,'PaperPositionMode','auto')
        print(fig,[figdir,'gradLDR_',datestr(startTime,'yyyymmdd_HHMMSS'),'_to_',datestr(endTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')
    end

    %% Plot powers and LDR

    disp('Plotting powers ...');

    fig=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',14);

    colormap('jet');

    s1=subplot(3,1,1);

    surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.DBMVC(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-105,-70]);
    colorbar

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('DBMVC (dB)')

    s2=subplot(3,1,2);
    colmapHX=jet(63);
    colmapHX=cat(1,[0.7,0,0],colmapHX);

    surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.DBMHX(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-106,-80]);
    colormap(s2,colmapHX);
    colorbar

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('DBMHX (dB)')

    s3=subplot(3,1,3);

    surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.LDR(:,plotInds),'EdgeColor','none');
    view(2)
    caxis([-25,0]);
    colorbar

    xlim([startTime,endTime]);
    ylim([0 ylimUpper]);

    ylabel('Altitude (km)')

    grid on
    box on

    title('LDR (dB)')

    linkaxes([s1 s2 s3],'xy')

    if saveFig
        set(gcf,'PaperPositionMode','auto')
        print(fig,[figdir,'powersLDR_',datestr(startTime,'yyyymmdd_HHMMSS'),'_to_',datestr(endTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')
    end
end