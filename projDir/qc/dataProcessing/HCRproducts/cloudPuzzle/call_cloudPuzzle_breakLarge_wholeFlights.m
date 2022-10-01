% Call cloud puzzle script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';
whichModel='era5';

plotFig=1;
showPlot='off';
saveData=0;
saveTime=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'cloudPuzzleEchoType/wholeFlights/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ=[];
    data.ECHO_TYPE_2D=[];
    data.FLAG=[];
    data.ANTFLAG=[];

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    ylimUpper=(max(data.asl(data.FLAG==1))./1000)+0.5;

    %% Truncate to non missing

    gapSecs=10;
    [data,nonMissingInds]=joinOverMissing(data,gapSecs);

    %% Create cloudID

    minCloudSizePix=1000;

    cloudID=makeCloudID(data,minCloudSizePix);

    %% Cloud puzzle
    % Breaks up really big clouds
    
    disp('Getting cloud puzzle ...')
    data.cloudPuzzle=f_cloudPuzzle_breakLarge(cloudID,data);

    uClouds=unique(data.cloudPuzzle(~isnan(data.cloudPuzzle)));
    cloudCount=length(uClouds);

    %% Un-truncate

    data=unJoinOverMissing(data,nonMissingInds);

    %% Plot
    if plotFig
        disp('Plotting ...');

        startPlot=startTime;

        while startPlot<endTime
            endPlot=startPlot+minutes(30);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);
            timeInds=timeInds(1:10:length(timeInds));

            timePlot=data.time(timeInds);
            dbzPlot=data.DBZ(:,timeInds);
            aslPlot=data.asl(:,timeInds);
            puzzlePlot=data.cloudPuzzle(:,timeInds);

            uClouds=unique(puzzlePlot(~isnan(puzzlePlot)));
            if isempty(uClouds)
                startPlot=endPlot;
                continue
            end

            cloudCount=length(uClouds);

            if cloudCount>0
                
                close all

                fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900],'visible',showPlot);

                %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ax1=subplot(2,1,1);
                hold on;

                sub1=surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
                view(2);
                sub1=colMapDBZ(sub1);
                ylim([-0.2,ylimUpper]);
                ylabel('Altitude (km)');
                xlim([timePlot(1),timePlot(end)]);
                title('Reflectivity')
                grid on

                ax2=subplot(2,1,2);

                colMapIn=jet(cloudCount);
                % Make order random
                indsCol=nan(cloudCount,1);
                indsCol(1:2:cloudCount)=1:2:cloudCount;
                reverse=fliplr(2:2:cloudCount);
                indsCol(2:2:cloudCount)=reverse;
                %indsCol=randperm(size(colMapIn,1));
                colMapInds=cat(2,indsCol,colMapIn);
                colMapInds=sortrows(colMapInds);
                colMap=colMapInds(:,2:end);
                
                hold on;
                sub2=surf(timePlot,aslPlot./1000,puzzlePlot,'edgecolor','none');
                view(2);
                ax2.Colormap=colMap;
                caxis([min(uClouds)-0.5 max(uClouds)+0.5])
                ylim([-0.2,ylimUpper]);
                ylabel('Altitude (km)');
                xlim([timePlot(1),timePlot(end)]);
                title('Cloud Puzzle')
                grid on
                cb=colorbar;
                cb.TickLabels={''};

                formatOut = 'yyyymmdd_HHMM';
                set(gcf,'PaperPositionMode','auto')
                print([figdir,project,'_Flight',num2str(aa),'_',datestr(timePlot(1),formatOut),'_to_',datestr(timePlot(end),formatOut),'_cloudPuzzle.png'],'-dpng','-r0');
            end
            startPlot=endPlot;
        end
    end
    %% Saving data

    if saveData
        disp('Saving cloud puzzle field ...')

        save([outdir,whichModel,'.cloudPuzzle.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'cloudPuzzle');

        if saveTime
            timeHCR=data.time;
            save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
        end
    end
    
end