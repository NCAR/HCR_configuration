% Compare data from before and after velocity correction

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
qcVersion='v1.2';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

saveTime=1;
plotYes=0;

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'velCorrPlots/wholeFlights/'];

for kk=1:size(caseList,1)

    disp(['Flight ',num2str(kk)]);
    disp('Loading HCR data.')

    startTime=datetime(caseList(kk,1:6));
    endTime=datetime(caseList(kk,7:12));

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    if ~isempty(fileList)
        data=[];

        data.VEL=[];

        dataVars=fieldnames(data);

        if length(fileList)==0
            disp('No data files found.');
            return
        end

        % Load data
        data=read_HCR(fileList,data,startTime,endTime);

        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end

        dataVars=dataVars(~cellfun('isempty',dataVars));

        %% Load from mat
        disp('Loading mat.');

        fileIn1=dir([modeldir,whichModel,'.velCorr.*.Flight',num2str(kk),'.mat']);
        velCorrIn=load([modeldir,fileIn1.name]);
        velCorrGet=velCorrIn.velCorr;

        fileIn2=dir([modeldir,whichModel,'.time.*.Flight',num2str(kk),'.mat']);
        velCorrTime=load([modeldir,fileIn2.name]);
        timeCorr=velCorrTime.timeHCR;

        %% Plot
        disp('Plotting ...');

        startPlot=startTime;

        while startPlot<endTime

            close all

            endPlot=startPlot+minutes(60);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);

            velMasked=data.VEL(:,timeInds);

            if sum(sum(~isnan(velMasked)))~=0

                timeIndsCorr=find(timeCorr>=data.time(timeInds(1)) & timeCorr<=data.time(timeInds(end)));

                velCorrPlot=velCorrGet(:,timeIndsCorr);
                
                %% Plot

                timeMasked=data.time(timeInds);
                aslMasked=data.asl(:,timeInds);
                elevMasked=data.elevation(:,timeInds);
                altMasked=data.altitude(:,timeInds);

                close all

                if etime(datevec(endPlot),datevec(startPlot))<=900
                    newInds=1:1:size(velMasked,2);
                elseif etime(datevec(endPlot),datevec(startPlot))<=3600
                    newInds=1:10:size(velMasked,2);
                else
                    newInds=1:100:size(velMasked,2);
                end

                % Resample for plotting
                newVEL=velMasked(:,newInds);
                newASL=aslMasked(:,newInds);
                newVelCorr=velCorrPlot(:,newInds);
                newTime=timeMasked(newInds);
                newElev=elevMasked(newInds);
                newAlt=altMasked(newInds);

                newVEL(:,newElev>0)=-newVEL(:,newElev>0);
                newVelCorr(:,newElev>0)=-newVelCorr(:,newElev>0);

                %% Plot hourly
                
                ylimupper=10;

                close all

                f1=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible','off');

                s1=subplot(2,1,1);
                hold on
                surf(newTime,newASL./1000,newVEL,'EdgeColor','none');
                view(2);
                plot(newTime,newAlt./1000,'-k','LineWidth',1.5);
                caxis([-5 5]);
                colM=colormap(velCols);
                colormap(s1,colM);
                colorbar
                title([project,' flight ',num2str(kk),' VEL'])
                ylim([0 ylimupper])
                xlim([newTime(1),newTime(end)])
                ylabel('Altitude (km)')
                grid on
                box on

                s2=subplot(2,1,2);
                hold on
                surf(newTime,newASL./1000,newVelCorr,'EdgeColor','none');
                view(2);
                plot(newTime,newAlt./1000,'-k','LineWidth',1.5);
                caxis([-5 5]);
                colM=colormap(velCols);
                colormap(s2,colM);
                colorbar
                title(['VEL CORR'])
                ylim([0 ylimupper])
                xlim([newTime(1),newTime(end)])
                ylabel('Altitude (km)')
                grid on
                box on

                set(gcf,'PaperPositionMode','auto')
                formatOut = 'yyyymmdd_HHMM';
                print(f1,[figdir,'velCorr_',datestr(newTime(1),formatOut),'_to_',datestr(newTime(end),formatOut)],'-dpng','-r0')
            end
            startPlot=endPlot;
        end
    end
end