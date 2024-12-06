% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';
whichModel='hrrr';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

saveData=1;
plotData=1;
ylimUpper=13;

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-14),'mergeLongShort/wholeIOPs/'];

snrThresh.DBZ=-6;
snrThresh.VEL=-5;
snrThresh.WIDTH=10;
snrThresh.LDRV=25;

%% Run processing

% Go through iops
for ii=2:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_short=[];
    data.VEL_unfold_short=[];
    data.WIDTH_short=[];
    data.SNRVC_short=[];
    data.LDRV_short=[];
    data.FLAG_short=[];

    data.DBZ_long=[];
    data.VEL_unfold_long=[];
    data.WIDTH_long=[];
    data.LDRV_long=[];
    data.FLAG_long=[];
       
    %% Load data
    disp('Loading data ...');

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
        if ~strcmp(thisName,'VEL')
            longField=data.([thisName,'_long']);
            shortField=data.([thisName,'_short']);
        else
            longField=data.VEL_unfold_long;
            shortField=data.VEL_unfold_short;
        end
        longField(data.FLAG_long~=1)=nan;
        shortField(data.FLAG_short~=1)=nan;

        % Mask short field with snr
        shortField(data.SNRVC_short<snrThresh.(thisName))=nan;

        % Fill in missing with long field
        outField=shortField;
        outField(isnan(shortField))=longField(isnan(shortField));

        % Output field
        if ~strcmp(thisName,'LDRV')
            dataOut.(thisName)=outField;
        else
            dataOut.LDR=outField;
        end
    end

     %% Save

    if saveData
        disp('Saving output fields.')

        dbz=dataOut.DBZ;
        save([outdir,whichModel,'.dbz.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(ii),'.mat'],'dbz','-v7.3');

        vel=dataOut.VEL;
        save([outdir,whichModel,'.vel.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(ii),'.mat'],'vel','-v7.3');

        width=dataOut.WIDTH;
        save([outdir,whichModel,'.width.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(ii),'.mat'],'width','-v7.3');

        ldr=dataOut.LDR;
        save([outdir,whichModel,'.ldr.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(ii),'.mat'],'ldr','-v7.3');
    end

    %% Plot in hourly increments

    if plotData
        disp('Plotting ...');

        startPlot=startTime;

        climsDbz=[-40,30];
        climsLdr=[-35,-5];
        climsVel=[-15,15];
        climsWidth=[0,3];

        while startPlot<endTime

            close all

            endPlot=startPlot+minutes(20);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);
            timeInds=timeInds(1:5:length(timeInds));

            timePlot=data.time(timeInds);
            dbzPlot=dataOut.DBZ(:,timeInds);
            velPlot=dataOut.VEL(:,timeInds);
            widthPlot=dataOut.WIDTH(:,timeInds);
            ldrPlot=dataOut.LDR(:,timeInds);

            if sum(sum(~isnan(dbzPlot)))~=0
                aslPlot=data.asl(:,timeInds);

                f1 = figure('Position',[200 500 1600 600],'DefaultAxesFontSize',12,'visible','off');
                t = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');

                s1=nexttile(1);
                surf(timePlot,aslPlot./1000,dbzPlot,'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                clim(climsDbz);
                s1.Colormap=jet;
                colorbar
                grid on
                box on
                title('Reflectivity (dBZ)');

                ylim([0 ylimUpper]);
                xlim([timePlot(1),timePlot(end)]);
                grid on
                box on

                s2=nexttile(2);
                surf(timePlot,aslPlot./1000,velPlot,'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                clim(climsVel);
                s2.Colormap=velCols;
                colorbar
                grid on
                box on
                title('Velocity (m s^{-1})');

                ylim([0 ylimUpper]);
                xlim([timePlot(1),timePlot(end)]);
                grid on
                box on

                s3=nexttile(3);
                surf(timePlot,aslPlot./1000,widthPlot,'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                clim(climsWidth);
                s3.Colormap=jet;
                colorbar
                grid on
                box on
                title('Spectrum width (m s^{-1})');

                ylim([0 ylimUpper]);
                xlim([timePlot(1),timePlot(end)]);
                grid on
                box on

                s4=nexttile(4);
                surf(timePlot,aslPlot./1000,ldrPlot,'edgecolor','none');
                view(2);
                ylabel('Altitude (km)');
                clim(climsLdr);
                s4.Colormap=jet;
                colorbar
                grid on
                box on
                title('Linear depol. ratio (dB)');

                ylim([0 ylimUpper]);
                xlim([timePlot(1),timePlot(end)]);
                grid on
                box on

                set(gcf,'PaperPositionMode','auto')
                print(f1,[figdir,project,'_IOP',num2str(ii),'_merged_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
            end
            startPlot=endPlot;
        end
    end
end
