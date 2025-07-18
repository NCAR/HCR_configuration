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

figdir=[indir(1:end-14),'mergeLongShort/calibrationBias/'];

%% Run processing

% Go through iops
for ii=1:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_short=[];
    data.VEL_short=[];
    data.WIDTH_short=[];
    data.SNRVC_short=[];
    data.LDRV_short=[];
    data.DBZ_long=[];
    data.VEL_long=[];
    data.WIDTH_long=[];
    data.LDRV_long=[];
       
    %% Load data
    disp('Loading data ...');

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.WIDTH_long(data.WIDTH_long<=0.1)=nan;
    data.WIDTH_short(data.WIDTH_short<=0.1)=nan;
    % data.DBZ_long(data.DBZ_long<=-40)=nan;
    % data.DBZ_short(data.DBZ_short<=-40)=nan;
    % data.LDRV_long(data.LDRV_long>-1)=nan;
    % data.LDRV_short(data.LDRV_short>-1)=nan;
    data.VEL_long(isnan(data.DBZ_long))=nan;
    data.VEL_short(isnan(data.DBZ_short))=nan;

    data.WIDTH_long(data.SNRVC_short<=10)=nan;
    data.WIDTH_short(data.SNRVC_short<=10)=nan;
    data.DBZ_long(data.SNRVC_short<=-6)=nan;
    data.DBZ_short(data.SNRVC_short<=-6)=nan;
    data.LDRV_long(data.SNRVC_short>25)=nan;
    data.LDRV_short(data.SNRVC_short>25)=nan;
    data.VEL_long(data.SNRVC_short<=-5)=nan;
    data.VEL_short(data.SNRVC_short<=-5)=nan;

    %% Sort by SNR
    disp('Loop through fields ...');

    dataFields={'DBZ','VEL','WIDTH','LDRV'};

    bySNR=[];
    for kk=1:size(dataFields,2)
        thisName=dataFields{kk};
        longField=data.([thisName,'_long']);
        shortField=data.([thisName,'_short']);
        shortLong=cat(2,shortField(:),longField(:));
        shortLong(any(isnan(shortLong),2),:)=[];
        bySNR.(thisName)=shortLong;
        if ii>1
            bySNRall.(thisName)=cat(1,bySNRall.(thisName),shortLong);
        end
    end

    if ii==1
        bySNRall=bySNR;
    end

    %% Plot

    edges.DBZ=[-100,-40:1:20,50];
    edges.VEL=[-30,-12:1:12,30];
    edges.WIDTH=[0:0.1:4,10];
    edges.LDRV=[-50,-30:1:10,30];

    f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12);
    t = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
    col=cat(1,[1,1,1],jet);
    
    for kk=1:size(dataFields,2)
         thisName=dataFields{kk};
       
        [N,~,~]=histcounts2(bySNR.(thisName)(:,1), ...
            bySNR.(thisName)(:,2),edges.(thisName),edges.(thisName));

        bias=mean(bySNR.(thisName)(:,1)-bySNR.(thisName)(:,2));
        rmse1=rmse(bySNR.(thisName)(:,1),bySNR.(thisName)(:,2));
        corr=corrcoef(bySNR.(thisName)(:,1),bySNR.(thisName)(:,2));

        plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

        s1=nexttile(kk);
        colormap(col);
        title([thisName,' (bias=',num2str(bias,2),', rmse=',num2str(rmse1,2),', corr=',num2str(corr(1,2),2),')'],'fontweight','bold','fontsize',12);
        hold on
        if sum(N,'all')>0
            pcolor(plotCoords,plotCoords,N);
            shading('flat');
            floorAx=prctile(bySNR.(thisName)(:),1);
            ceilAx=prctile(bySNR.(thisName)(:),99);
            xlim([floorAx,ceilAx]);
            ylim([floorAx,ceilAx]);
            plot([floorAx,ceilAx],[floorAx,ceilAx],'-m','LineWidth',2);
            box on
            xlabel('Long pulse')
            ylabel('Short pulse')
            s1.SortMethod='childorder';
            set(gca,'layer','top')
        end
        
    end
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_IOP',num2str(ii)],'-dpng','-r0')
end

%% Plot all iops together

% Remove aliased
bySNRall.VELorig=bySNRall.VEL;
bySNRall.VEL(abs(bySNRall.VEL(:,1)-bySNRall.VEL(:,2))>7.8,:)=[];

f1 = figure('Position',[200 500 800 700],'DefaultAxesFontSize',12);
t = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
col=cat(1,[1,1,1],jet);

for kk=1:size(dataFields,2)
    thisName=dataFields{kk};

    [N,~,~]=histcounts2(bySNRall.(thisName)(:,1), ...
        bySNRall.(thisName)(:,2),edges.(thisName),edges.(thisName));

    bias=mean(bySNRall.(thisName)(:,1)-bySNRall.(thisName)(:,2));
    rmse1=rmse(bySNRall.(thisName)(:,1),bySNRall.(thisName)(:,2));
    corr=corrcoef(bySNRall.(thisName)(:,1),bySNRall.(thisName)(:,2));

    plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

    s1=nexttile(kk);
    colormap(col);
    title([thisName,' (bias=',num2str(bias,2),', rmse=',num2str(rmse1,2),', corr=',num2str(corr(1,2),2),')'],'fontweight','bold','fontsize',12);
    hold on
    if sum(N,'all')>0
        pcolor(plotCoords,plotCoords,N);
        shading('flat');
        floorAx=prctile(bySNRall.(thisName)(:),0.1);
        ceilAx=prctile(bySNRall.(thisName)(:),99.9);
        xlim([floorAx,ceilAx]);
        ylim([floorAx,ceilAx]);
        plot([floorAx,ceilAx],[floorAx,ceilAx],'-m','LineWidth',2);
        box on
        xlabel('Long pulse')
        ylabel('Short pulse')
        s1.SortMethod='childorder';
        set(gca,'layer','top')
    end

end
set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project],'-dpng','-r0')