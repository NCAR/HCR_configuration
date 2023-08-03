% Compare data from before and after velocity correction

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

saveTime=0;
plotYes=1;
showPlot='off';

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'velCorrZenithPlots/wholeFlights/'];

if strcmp(project,'socrates')
    polyTimePeriod=20; %Time period for poly fit in seconds
else
    polyTimePeriod=15; %Time period for poly fit in seconds
end
polyOrder=3; % Order of polynomial fit

for kk=1:size(caseList,1)
    
    disp(['Flight ',num2str(kk)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(kk,1:6));
    endTime=datetime(caseList(kk,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        data=[];
        
        data.vertical_velocity=[];
        data.nyquist_velocity=[];
        data.DBZ=[];
        data.VEL=[];
        data.TOPO=[];
        data.FLAG=[];
        data.ANTFLAG=[];
        
        dataVars=fieldnames(data);
              
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
                    
        %% Correct nadir pointing vel
        disp('Nadir correction')
        [nadirCorrection,surfInd]=velCorrNadir(data,polyTimePeriod,polyOrder);

        % Vel corr nadir
        velCorrN=data.VEL-nadirCorrection';
        velNadirCorrSurf=velCorrN(surfInd);

        % Corrected nadir vel
        velNadirCorr=data.VEL;
        velNadirCorr(:,data.elevation<=0)=velCorrN(:,data.elevation<=0);

        %% Correct zenith pointing vel
        disp('Zenith correction')
        %[zenithCorrection,medCloudVelNadir,medCloudVelZenith,smoothfactor,velZenithCorrSmooth]=velCorrZenith(data,velNadirCorr);
        [zenithCorrection,smoothfactor]=velCorrZenith(data,velNadirCorr);

        % Vel corr zenith
        velCorrZ=data.VEL-zenithCorrection;

        % Corrected nadir vel
        velAllCorr=velNadirCorr;
        velAllCorr(:,data.elevation>0)=velCorrZ(:,data.elevation>0);

        %% Save
        disp('Saving velCorr field.')

        velCorr=velAllCorr;
        save([outdir,whichModel,'.velCorr.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(kk),'.mat'],'velCorr');

        if saveTime
            timeHCR=data.time;
            save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
                datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(kk),'.mat'],'timeHCR');
        end

        %% Plot

        if plotYes

            disp('Plotting.');

            plotInds=1:10:length(data.time);
            plotVelMasked=data.VEL;
            plotVelMasked(data.FLAG~=1)=nan;
            plotVelMasked(:,data.ANTFLAG>2)=nan;
            plotVelMasked(:,data.elevation>0)=-plotVelMasked(:,data.elevation>0);

            plotVel=plotVelMasked(:,plotInds);
            plotAsl=data.asl(:,plotInds);
            plotTime=data.time(plotInds);

            velAllPlot=velAllCorr;
            velAllPlot(data.FLAG~=1)=nan;
            velAllPlot(:,data.ANTFLAG>2)=nan;
            velAllPlot(:,data.elevation>0)=-velAllPlot(:,data.elevation>0);
            plotVelCorr=velAllPlot(:,plotInds);

            aslTest=data.asl;
            aslTest(isnan(plotVelMasked))=nan;
            ylimupper=ceil(max(max(aslTest,[],'omitnan')))./1000;

            velNadirCorrSmooth=movmedian(velNadirCorrSurf,smoothfactor,'omitnan');
            velNadirCorrSmooth(isnan(velNadirCorrSurf))=nan;
            zenithCorrection(data.elevation<=0)=nan;

            close all

            f1=figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible',showPlot);

            s1=subplot(3,1,1);
            hold on
            surf(plotTime,plotAsl./1000,plotVel,'EdgeColor','none');
            view(2);
            plot(data.time,data.altitude./1000,'-k','LineWidth',1.5);
            caxis([-5 5]);
            colM=colormap(velCols);
            colormap(s1,colM);
            colorbar
            title([project,' flight ',num2str(kk),' VEL'])
            ylim([0 ylimupper])
            xlim([data.time(1),data.time(end)])
            ylabel('Altitude (km)')
            grid on
            box on

            s2=subplot(3,1,2);
            hold on
            surf(plotTime,plotAsl./1000,plotVelCorr,'EdgeColor','none');
            view(2);
            plot(data.time,data.altitude./1000,'-k','LineWidth',1.5);
            caxis([-5 5]);
            colM=colormap(velCols);
            colormap(s2,colM);
            colorbar
            title(['VEL CORR'])
            ylim([0 ylimupper])
            xlim([data.time(1),data.time(end)])
            ylabel('Altitude (km)')
            grid on
            box on

            s3=subplot(3,1,3);
            %plot(data.time,meanCloudVel,'-b','LineWidth',0.8);
            hold on
            plot(data.time,velNadirCorrSmooth,'-y','LineWidth',1.5);
            plot(data.time,nadirCorrection,'-g','LineWidth',1);
            plot(data.time,zenithCorrection,'-c','LineWidth',1.5);
            ylim([-1.5 1.5])
            xlim([data.time(1),data.time(end)])
            ylabel('Velocities (m/s)')
            legend('Surface vel','Nadir correction','Zenith correction','Orientation','horizontal')
            grid on
            box on
            s3.SortMethod='childorder';

            s1pos=s1.Position;
            s2pos=s2.Position;
            s3pos=s3.Position;
            s1.Position=[s1pos(1),s1pos(2),s3pos(3),s1pos(4)];
            s2.Position=[s2pos(1),s2pos(2),s3pos(3),s2pos(4)];

            set(gcf,'PaperPositionMode','auto')
            print(f1,[figdir,project,'_rf',num2str(kk),'cloudTop_vel'],'-dpng','-r0')
        end
    end
end