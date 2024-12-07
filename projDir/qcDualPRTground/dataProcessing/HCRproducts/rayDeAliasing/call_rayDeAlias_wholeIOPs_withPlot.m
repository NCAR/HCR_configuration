% Plot HCR pid from mat file in hourly plots

% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz_combined'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';
whichModel='hrrr';
dataFreq=10;

saveData=1;
plotData=1;

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

ylimUpper=13;

nyquistFile='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/moments/50hz_longPulse/20240510/cfrad.20240510_180000.031_to_20240510_180100.028_HCR.nc';
readNV=ncread(nyquistFile,'nyquist_velocity');
nyq=mode(readNV);

indir=HCRdir(project,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

figdir=[indir(1:end-14),'rayDeAlias/wholeIOPs/'];

for aa=1:size(caseList,1)

    disp(['IOP ',num2str(aa)]);
    disp('Loading HCR data.')
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    data=[];
    data.VEL_long=[];
    data.FLAG_long=[];
    data.VEL_unfold_short=[];
    data.FLAG_short=[];
                
    dataVars=fieldnames(data);
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
       
    velMaskedLong=data.VEL_long;
    velMaskedLong(data.FLAG_long~=1)=nan;

    velMaskedUnfold=data.VEL_unfold_short;
    velMaskedUnfold(data.FLAG_short~=1)=nan;

    %% De-alias with unfolded

    checkFold=[2,4,6];

    velDeAliasCheck=velMaskedUnfold-velMaskedLong;

    deAliasMaskE=zeros(size(velDeAliasCheck));
    for jj=1:3
        deAliasMaskE(velDeAliasCheck>checkFold(jj)*nyq-3)=checkFold(jj)*nyq;
        deAliasMaskE(velDeAliasCheck<-(checkFold(jj)*nyq-3))=-checkFold(jj)*nyq;
    end

    velLongTemp=velMaskedLong+deAliasMaskE;
   
    %% Correct velocity folding
  
    disp('De-aliasing ...');
    velDeAliased=dealiasByRay(velLongTemp,data.elevation,nyq,dataFreq,data.time,[]);

    %% De-alias with unfolded again

    velDeAliasCheck=velMaskedUnfold-velDeAliased;

    deAliasMaskE=zeros(size(velDeAliasCheck));
    for jj=1:3
        deAliasMaskE(velDeAliasCheck>checkFold(jj)*nyq-3)=checkFold(jj)*nyq;
        deAliasMaskE(velDeAliasCheck<-(checkFold(jj)*nyq-3))=-checkFold(jj)*nyq;
    end

    % Remove small regions
    uMask=unique(deAliasMaskE);
    for kk=1:length(uMask)
        if uMask(kk)~=0
            thisMask=deAliasMaskE==uMask(kk);
            thisMaskShrink=bwareaopen(thisMask,1000);
            deAliasMaskE(thisMaskShrink==0 & thisMask==1)=0;
        end
    end

    velDeAliased=velDeAliased+deAliasMaskE;
   
    %% Save

    if saveData
        disp('Saving velUnfolded field.')

        velUnfolded=velDeAliased;

        save([outdir,whichModel,'.velUnfolded_long.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'velUnfolded','-v7.3');
    end

    %% Plot in hourly increments

    if plotData
        disp('Plotting ...');

        startPlot=startTime;

        while startPlot<endTime

            close all

            endPlot=startPlot+minutes(20);
            timeInds=find(data.time>=startPlot & data.time<=endPlot);
            timeInds=timeInds(1:3:length(timeInds));

            timePlot=data.time(timeInds);
            velPlot=data.VEL_long(:,timeInds);

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
                    clim([-13 13]);
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
                    clim([-13 13]);
                    ylim([0 ylimUpper]);
                    xlim([timePlot(1),timePlot(end)]);
                    colorbar
                    grid on
                    box on

                    set(gcf,'PaperPositionMode','auto')
                    print(f1,[figdir,project,'_IOP',num2str(aa),'_velFinal_',datestr(timePlot(1),'yyyymmdd_HHMMSS'),'_to_',datestr(timePlot(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
                end
            end
            startPlot=endPlot;
        end
    end
end