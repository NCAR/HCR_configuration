% Plot HCR pid from mat file in hourly plots

% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
dataFreq=10;
whichModel='era5';

saveData=1;
plotData=1;

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

ylimUpper=13;

indir=HCRdir(project,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

figdir=[indir(1:end-5),'rayDeAlias/wholeFlights/'];

for aa=1:size(caseList,1)

    disp(['Flight ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    data.nyquist_velocity=[];
    data.VEL_CORR=[];
    data.FLAG=[];
                
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
       
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end

    dataVars=dataVars(~cellfun('isempty',dataVars));

    velMasked=data.VEL_CORR;
    velMasked(data.FLAG~=1)=nan;

    %% Correct velocity folding

    nyq=mode(data.nyquist_velocity);
  
    disp('De-aliasing ...');
    velDeAliased=dealiasByRay(velMasked,data.elevation,nyq,dataFreq,data.time,[]);

    %% Save

    if saveData
        disp('Saving velUnfolded field.')

        velUnfolded=velDeAliased;

        save([outdir,whichModel,'.velUnfolded.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'velUnfolded');
    end

    %% Plot in hourly increments

    if plotData
        disp('Plotting ...');

        startPlot=startTime;

        while startPlot<endTime

            close all

            endPlot=startPlot+minutes(15);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);

            timePlot=data.time(timeInds);
            velPlot=data.VEL_CORR(:,timeInds);

            if sum(sum(~isnan(velPlot)))~=0
                aslPlot=data.asl(:,timeInds);
                velFinalPlot=velDeAliased(:,timeInds);

                if sum(sum(~isnan(velFinalPlot)))~=0

                    f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'visible','off');

                    s1=subplot(2,1,1);

                    colormap(velCols)

                    hold on
                    surf(timePlot,aslPlot./1000,velPlot,'edgecolor','none');
                    view(2);
                    ylabel('Altitude (km)');
                    caxis([-13 13]);
                    ylim([0 ylimUpper]);
                    xlim([timePlot(1),timePlot(end)]);
                    colorbar
                    grid on
                    box on
                    title('VEL CORR')

                    s2=subplot(2,1,2);

                    hold on
                    surf(timePlot,aslPlot./1000,velFinalPlot,'edgecolor','none');
                    view(2);
                    ylabel('Altitude (km)');
                    title(['VEL UNFOLDED']);
                    caxis([-13 13]);
                    ylim([0 ylimUpper]);
                    xlim([timePlot(1),timePlot(end)]);
                    colorbar
                    grid on
                    box on

                    set(gcf,'PaperPositionMode','auto')
                    print(f1,[figdir,project,'_Flight',num2str(aa),'_velFinal_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
                end
            end
            startPlot=endPlot;
        end
    end
end