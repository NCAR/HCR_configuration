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

figdir=[indir(1:end-14),'mergeLongShort/analyse/'];

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

    %% Mask short
    dataFields={'DBZ','VEL','WIDTH','LDRV'};

    bySNR=[];
    for jj=1:length(snrEdges.DBZ)-1
        for kk=1:size(dataFields,2)
            thisName=dataFields{kk};
            thisMask=(data.SNRVC_short>snrEdges.(thisName)(jj) & data.SNRVC_short<=snrEdges.(thisName)(jj+1));
            longField=data.([thisName,'_long']);
            shortField=data.([thisName,'_short']);
            shortLong=cat(2,shortField(thisMask==1),longField(thisMask==1));
            shortLong(any(isnan(shortLong),2),:)=[];
            bySNR.(['bin',num2str(jj)]).(thisName)=shortLong;
            if ii>1
                bySNRall.(['bin',num2str(jj)]).(thisName)=cat(1,bySNRall.(['bin',num2str(jj)]).(thisName),shortLong);
            end
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

    for kk=1:size(dataFields,2)
        close all
        f1 = figure('Position',[200 500 1600 600],'DefaultAxesFontSize',12);
        t = tiledlayout(2,5,'TileSpacing','tight','Padding','compact');
        col=cat(1,[1,1,1],jet);
        colormap(col);

        thisName=dataFields{kk};
        title(t,[thisName],'fontweight','bold','fontsize',16);
        for jj=1:length(snrEdges.DBZ)-1
            [N,~,~]=histcounts2(bySNR.(['bin',num2str(jj)]).(thisName)(:,1), ...
                bySNR.(['bin',num2str(jj)]).(thisName)(:,2),edges.(thisName),edges.(thisName));

            plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

            s1=nexttile(jj);
            hold on
            if sum(N,'all')>0
                pcolor(plotCoords,plotCoords,N);
                shading('flat');
                floorAx=prctile(bySNR.(['bin',num2str(jj)]).(thisName)(:),1);
                ceilAx=prctile(bySNR.(['bin',num2str(jj)]).(thisName)(:),99);
                xlim([floorAx,ceilAx]);
                ylim([floorAx,ceilAx]);
                plot([floorAx,ceilAx],[floorAx,ceilAx],'-m','LineWidth',2);
                box on
                xlabel('Long pulse')
                ylabel('Short pulse')
                s1.SortMethod='childorder';
                set(gca,'layer','top')
            end
            title([num2str(snrEdges.(thisName)(jj)),' to ',num2str(snrEdges.(thisName)(jj+1)),' dB SNR']);
        end
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_IOP',num2str(ii),'_',thisName],'-dpng','-r0')
    end
end

%% Plot all iops together
for kk=1:size(dataFields,2)
    close all
    f1 = figure('Position',[200 500 1600 600],'DefaultAxesFontSize',12);
    t = tiledlayout(2,5,'TileSpacing','tight','Padding','compact');
    col=cat(1,[1,1,1],jet);
    colormap(col);

    thisName=dataFields{kk};
    title(t,[thisName],'fontweight','bold','fontsize',16);
    for jj=1:length(snrEdges.DBZ)-1
        [N,~,~]=histcounts2(bySNRall.(['bin',num2str(jj)]).(thisName)(:,1), ...
            bySNRall.(['bin',num2str(jj)]).(thisName)(:,2),edges.(thisName),edges.(thisName));

        plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

        s1=nexttile(jj);
        hold on
        if sum(N,'all')>0
            pcolor(plotCoords,plotCoords,N);
            shading('flat');
            floorAx=prctile(bySNRall.(['bin',num2str(jj)]).(thisName)(:),1);
            ceilAx=prctile(bySNRall.(['bin',num2str(jj)]).(thisName)(:),99);
            xlim([floorAx,ceilAx]);
            ylim([floorAx,ceilAx]);
            plot([floorAx,ceilAx],[floorAx,ceilAx],'-m','LineWidth',2);
            box on
            xlabel('Long pulse')
            ylabel('Short pulse')
            s1.SortMethod='childorder';
            set(gca,'layer','top')
        end
        title([num2str(snrEdges.(thisName)(jj)),' to ',num2str(snrEdges.(thisName)(jj+1)),' dB SNR']);
    end
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_',thisName],'-dpng','-r0')
end