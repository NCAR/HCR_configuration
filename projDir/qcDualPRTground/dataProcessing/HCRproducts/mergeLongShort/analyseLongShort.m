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
    data.NCP_short=[];
    data.LDRV_short=[];
    data.DBZ_long=[];
    data.VEL_long=[];
    data.WIDTH_long=[];
    data.SNRVC_long=[];
    data.NCP_long=[];
    data.LDRV_long=[];
       
    %% Load data
    disp('Loading data ...');

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.WIDTH_long(data.WIDTH_long<=0.1)=nan;
    data.WIDTH_short(data.WIDTH_short<=0.1)=nan;
    data.DBZ_long(data.DBZ_long<=-40)=nan;
    data.DBZ_short(data.DBZ_short<=-40)=nan;

    %% Sort by SNR
    disp('Sorting by SNR ...');
    snrEdges=[-inf,-20,-15,-10:5,10:10:50,inf];

    dataFields={'DBZ','VEL','WIDTH','SNRVC','NCP','LDRV'};

    bySNR=[];
    for jj=1:length(snrEdges)-1
        thisMask=(data.SNRVC_short>snrEdges(jj) & data.SNRVC_short<=snrEdges(jj+1));
        for kk=1:size(dataFields,2)
            thisName=dataFields{kk};
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
    edges.SNRVC=[-50,-20:5:60,100];
    edges.NCP=[-1,0:0.05:1,2];
    edges.LDRV=[-50,-30:1:10,30];

    for jj=1:length(snrEdges)-1
        close all
        f1 = figure('Position',[200 500 1200 750],'DefaultAxesFontSize',12);
        t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
        col=cat(1,[1,1,1],jet);
        colormap(col);
        title(t,[num2str(snrEdges(jj)),' to ',num2str(snrEdges(jj+1)),' dB SNR'],'fontweight','bold','fontsize',16);
        for kk=1:size(dataFields,2)
            thisName=dataFields{kk};
            [N,~,~]=histcounts2(bySNR.(['bin',num2str(jj)]).(thisName)(:,1), ...
                bySNR.(['bin',num2str(jj)]).(thisName)(:,2),edges.(thisName),edges.(thisName));

            plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

            s1=nexttile(kk);
            if sum(N,'all')>0
                pcolor(plotCoords,plotCoords,N);
                shading('flat');
                colorbar
                floorAx=prctile(bySNR.(['bin',num2str(jj)]).(thisName)(:),1);
                ceilAx=prctile(bySNR.(['bin',num2str(jj)]).(thisName)(:),99);
                xlim([floorAx,ceilAx]);
                ylim([floorAx,ceilAx]);
                box on
                xlabel('Long pulse')
                ylabel('Short pulse')
                set(gca,'layer','top')
            end
            title([thisName,', ',num2str(sum(N,'all')),' samples']);
        end
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_IOP',num2str(ii),'_bin',num2str(jj)],'-dpng','-r0')
    end
end

%% Plot all iops together
for jj=1:length(snrEdges)-1
    close all
    f1 = figure('Position',[200 500 1200 750],'DefaultAxesFontSize',12);
    t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
    col=cat(1,[1,1,1],jet);
    colormap(col);
    title(t,[num2str(snrEdges(jj)),' to ',num2str(snrEdges(jj+1)),' dB SNR'],'fontweight','bold','fontsize',16);
    for kk=1:size(dataFields,2)
        thisName=dataFields{kk};
        [N,~,~]=histcounts2(bySNRall.(['bin',num2str(jj)]).(thisName)(:,1), ...
            bySNRall.(['bin',num2str(jj)]).(thisName)(:,2),edges.(thisName),edges.(thisName));

        plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;

        s1=nexttile(kk);
        if sum(N,'all')>0
            pcolor(plotCoords,plotCoords,N);
            shading('flat');
            colorbar
            floorAx=prctile(bySNRall.(['bin',num2str(jj)]).(thisName)(:),1);
            ceilAx=prctile(bySNRall.(['bin',num2str(jj)]).(thisName)(:),99);
            xlim([floorAx,ceilAx]);
            ylim([floorAx,ceilAx]);
            box on
            xlabel('Long pulse')
            ylabel('Short pulse')
            set(gca,'layer','top')
        end
        title([thisName,', ',num2str(sum(N,'all')),' samples']);
    end
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_bin',num2str(jj)],'-dpng','-r0')
end