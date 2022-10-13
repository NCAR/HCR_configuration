% Call cloud classification script

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='noreaster'; %socrates, aristo, cset, otrec
quality='qc2'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v2.0';
whichModel='era5';

plotFig=1;
showPlot='off';
saveCloudClass=1;
saveCloudPuzzle=1;
saveTime=0;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'cloudClass/wholeFlights/'];

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
    
    data.DBZ=[];
    data.ECHO_TYPE_2D=[];
    data.FLAG=[];
    data.ANTFLAG=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    data.TOPO=[];
        
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
 
    ylimUpper=(max(data.asl(~isnan(data.ECHO_TYPE_2D)))./1000)+0.5;
   
     %% Truncate to non missing

    gapSecs=10;
    [data,nonMissingInds]=joinOverMissing(data,gapSecs);

    %% Create cloudID

    minCloudSizePix=1000;

    cloudID=makeCloudID(data,minCloudSizePix);

    %% Cloud puzzle
    % Breaks up really big clouds
    
    disp('Breaking up large ...')
    data.cloudPuzzle=f_cloudPuzzle_breakLarge(cloudID,data);

    % Breaks out isolated convective that penetrate into stratiform
    disp('Breaking out isolated conv ...')
    data.cloudPuzzle_echoType=breakout_isolatedConv(data);

    %% Un-truncate

    disp('Adding gaps back in ...')
    data=unJoinOverMissing(data,nonMissingInds);
       
    %% Cloud classification

    disp('Cloud classification ...')
    cloudClass=findCloudClass(data.ECHO_TYPE_2D,data.cloudPuzzle_echoType,data.TEMP,data.MELTING_LAYER,data.elevation,data.TOPO,data.asl);

    % Fill in small with not classified
    cloudClass(~isnan(data.ECHO_TYPE_2D) & isnan(cloudClass))=0;

    %% Plot
    if plotFig

        colmapCC=[0,0,0;
            204,255,204;
            153,204,0;
            0,128,0;
            0,204,255;
            51,102,255;
            0,0,180;
            255,204,0;
            255,102,0;
            220,0,0;
            255,153,220;
            204,153,255;
            128,0,128];

        colmapCC=colmapCC./255;

        colmapSC=[0,0.1,0.6;
            0.38,0.42,0.96;
            0.65,0.74,0.86;
            0.32,0.78,0.59;
            0.7,0,0;
            1,0,1;
            1,1,0;
            0.99,0.77,0.22;
            1,0,0];

        disp('Plotting ...');

        startPlot=startTime;

        while startPlot<endTime
            endPlot=startPlot+minutes(30);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);
            timeInds=timeInds(1:10:length(timeInds));

            plotTime=data.time(timeInds);
            plotAsl=data.asl(:,timeInds);

            classPlot=cloudClass(:,timeInds);
            classPlot(classPlot==11)=4;
            classPlot(classPlot==12)=5;
            classPlot(classPlot==13)=6;
            classPlot(classPlot==21)=7;
            classPlot(classPlot==22)=8;
            classPlot(classPlot==23)=9;
            classPlot(classPlot==31)=10;
            classPlot(classPlot==32)=11;
            classPlot(classPlot==33)=12;

            if sum(sum(~isnan(classPlot)))==0
                startPlot=endPlot;
                continue
            end

            csSubPlot=data.ECHO_TYPE_2D(:,timeInds);
            csSubPlot(csSubPlot==14)=1;
            csSubPlot(csSubPlot==16)=2;
            csSubPlot(csSubPlot==18)=3;
            csSubPlot(csSubPlot==25)=4;
            csSubPlot(csSubPlot==30)=5;
            csSubPlot(csSubPlot==32)=6;
            csSubPlot(csSubPlot==34)=7;
            csSubPlot(csSubPlot==36)=8;
            csSubPlot(csSubPlot==38)=9;

            close all

            f1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900],'visible',showPlot);
            s1=subplot(2,1,1);

            hold on
            surf(plotTime,plotAsl./1000,csSubPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([0 10]);
            ylim([0 ylimUpper]);
            xlim([plotTime(1),plotTime(end)]);
            s1.Colormap=colmapSC;
            caxis([0.5 9.5]);
            cb=colorbar;
            cb.Ticks=1:9;
            cb.TickLabels={'Strat Low','Strat Mid','Strat High','Mixed',...
                'Conv','Conv Elev','Conv Shallow','Conv Mid','Conv Deep'};
            grid on
            box on
            title('Echo type')

            s2=subplot(2,1,2);

            hold on
            surf(plotTime,plotAsl./1000,classPlot,'edgecolor','none');
            view(2);
            ylabel('Altitude (km)');
            caxis([-1 13]);
            ylim([0 ylimUpper]);
            xlim([plotTime(1),plotTime(end)]);
            s2.Colormap=colmapCC;
            caxis([-0.5 12.5]);
            cb=colorbar;
            cb.Ticks=0:12;
            cb.TickLabels={'Not classified','Strat Low','Strat Mid','Strat High',...
                'Strat Precip Shallow','Strat Precip Mid','Strat Precip Deep',...
                'Conv Young Shallow','Conv Young Mid','Conv Yong Deep',...
                'Conv Mature Shallow','Conv Mature Mid','Conv Mature Deep'};
            grid on
            box on
            title('Cloud classification')

            linkaxes([s1],'xy');

            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_cloudClass_',datestr(plotTime(1),'yyyymmdd_HHMMSS'),'_to_',datestr(plotTime(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
            startPlot=endPlot;
        end
    end

    %% Saving data

    cloudPuzzle=data.cloudPuzzle_echoType;

    if saveCloudClass

        disp('Saving cloudClass field ...')

        save([outdir,whichModel,'.cloudClass.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'cloudClass');
    end

    if saveCloudPuzzle

        disp('Saving cloudPuzzle field ...')

        save([outdir,whichModel,'.cloudPuzzle.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'cloudPuzzle');
    end

    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end

end