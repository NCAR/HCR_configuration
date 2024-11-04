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

snrThresh.DBZ=-6;
snrThresh.VEL=-5;
snrThresh.WIDTH=10;
snrThresh.LDRV=25;

%% Run processing

% Go through iops
for ii=1:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_short=[];
    data.VEL_unfold=[];
    data.WIDTH_short=[];
    data.SNRVC_short=[];
    data.LDRV_short=[];
    data.FLAG_short=[];

    data.DBZ_long=[];
    data.VEL_long=[];
    data.WIDTH_long=[];
    data.LDRV_long=[];
    data.FLAG_long=[];
       
    %% Load data
    disp('Loading data.FLAG_short=[];data ...');

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    %% Loop through fields
    dataFields={'DBZ','VEL','WIDTH','LDRV'};

    dataOut=[];
    for kk=1:size(dataFields,2)
        thisName=dataFields{kk};

        % Mask fields with flag fields
        longField=data.([thisName,'_long']);
        longField(data.FLAG_long~=1)=nan;
        if ~strcmp(thisName,'VEL')
            shortField=data.([thisName,'_short']);
        else
            shortField=data.VEL_unfold;
        end
        shortField(data.FLAG_short~=1)=nan;

        % Mask short field with snr
        if ~strcmp(thisName,'VEL')
            shortField(data.SNRVC<snrThresh.(thisName))=nan;
        end

        % Output field
        dataOut.(thisName)=shortField;

        % Fill in missing with long field
        dataOut.(thisName)(isnan(shortField))=longField(isnan(shortField));
        
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