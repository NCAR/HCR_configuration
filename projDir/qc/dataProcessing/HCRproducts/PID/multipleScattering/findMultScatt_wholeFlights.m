% Plot HCR variables

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

indir=HCRdir(project,quality,qcVersion,freqData);

ylimUpper=15;

saveFig=1;

figdir=[indir(1:end-5),'multScatt/wholeFlights/'];
if ~exist(figdir, 'dir')
    mkdir(figdir)
end

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    disp('Loading data ...');

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ_MASKED=[];
    data.DBMVC=[];
    data.DBMHX=[];
    data.LDR=[];
    data.ICING_LEVEL=[];

    data=read_HCR(fileList,data,startTime,endTime);

    data.DBMHX(isnan(data.LDR))=-106;
    data.DBMHX(isnan(data.DBZ_MASKED))=nan;

    data.DBMVC(isnan(data.DBZ_MASKED))=nan;
    data.LDR(isnan(data.DBZ_MASKED))=nan;

    %% Find icing regions
    [LDRgrad,smoothLDR,iceMask]=findIce(data);

    %% Plot

    startPlot=startTime;

    if saveFig
        disp('Plotting ...');
        while startPlot<endTime
            endPlot=startPlot+minutes(30);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);
            plotInds=timeInds(1:10:length(timeInds));

            if sum(sum(~isnan(data.LDR(:,plotInds))))==0
                startPlot=endPlot;
                continue
            end

            %% Plot LDR, smooth and gradient
            close all
            
            fig=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',14,'visible','off');

            colormap('jet');

            s1=subplot(3,1,1);

            surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.LDR(:,plotInds),'EdgeColor','none');
            view(2)
            caxis([-25,0]);
            colorbar

            xlim([startPlot,endPlot]);
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

            xlim([startPlot,endPlot]);
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

            iceMaskPlot=iceMask(:,plotInds);
            timePlot=data.time(plotInds);
            aslPlot=data.asl(:,plotInds);
            bounds=bwboundaries(iceMaskPlot);
            for kk=1:length(bounds)
                boundary=bounds{kk};
                btimes=timePlot(boundary(:,2));
                linInds=sub2ind(size(aslPlot),boundary(:,1),boundary(:,2));
                basl=aslPlot(linInds)./1000;
                plot(btimes,basl,'k','LineWidth',1.5);
            end

            s3.SortMethod='childorder';

            xlim([startPlot,endPlot]);
            ylim([0 ylimUpper]);

            ylabel('Altitude (km)')

            grid on
            box on

            title('LDR gradient (dB)')

            linkaxes([s1 s2 s3],'xy')

            set(gcf,'PaperPositionMode','auto')
            print(fig,[figdir,'gradLDR_Flight',num2str(aa),'_',datestr(startPlot,'yyyymmdd_HHMMSS'),'_to_',datestr(endPlot,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')

            %% Plot powers and LDR

            fig=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',14,'visible','off');

            colormap('jet');

            s1=subplot(3,1,1);

            surf(data.time(plotInds),data.asl(:,plotInds)./1000,data.DBMVC(:,plotInds),'EdgeColor','none');
            view(2)
            caxis([-105,-70]);
            colorbar

            xlim([startPlot,endPlot]);
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

            xlim([startPlot,endPlot]);
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

            xlim([startPlot,endPlot]);
            ylim([0 ylimUpper]);

            ylabel('Altitude (km)')

            grid on
            box on

            title('LDR (dB)')

            linkaxes([s1 s2 s3],'xy')

            set(gcf,'PaperPositionMode','auto')
            print(fig,[figdir,'powersLDR_Flight',num2str(aa),'_',datestr(startPlot,'yyyymmdd_HHMMSS'),'_to_',datestr(endPlot,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0')

            startPlot=endPlot;
        end
    end
end