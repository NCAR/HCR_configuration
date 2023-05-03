% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

showPlot='off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

if strcmp(project,'cset') | strcmp(project,'socrates') | ...
        strcmp(project,'otrec') | strcmp(project,'spicule')
    upperLimDBZ=12;
    upperLimVEL=5;
elseif strcmp(project,'noreaster')
    upperLimDBZ=14;
    upperLimVEL=7;
else
    error('Set upperLimDBZ and upperLimVE')
end

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'attenPlots/wholeFlights/'];

%[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp(['Starting at ',datestr(datetime('now'),'yyyy-mm-dd HH:MM')]);
    disp('Loading data ...');

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get data

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ = [];
    data.U_SURF=[];
    data.V_SURF=[];
    data.SST=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    data.ANTFLAG=[];
    data.rotation=[];
    %data.MELTING_LAYER=[];
    data.pulse_width=[];

    dataVars=fieldnames(data);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.frq=ncread(fileList{1},'frequency');

    %% Correct for gaseous attenuation

    disp('Calculating gaseous attenuation ...');
    [~,gasAttCloud,~,gasAttCloudMat]=get_gas_atten(data);

    % Extend to surface
    [linInd,~,~]=hcrSurfInds(data);
    for jj=1:2
        replaceLinInd=linInd;
        replaceLinInd(isnan(linInd) | ~isnan(gasAttCloudMat(linInd)) | data.FLAG(linInd)~=7)=[];

        gasAttCloudMat(replaceLinInd)=gasAttCloudMat(replaceLinInd-jj);
        if jj==2
            gasAttCloudMat(replaceLinInd-1)=gasAttCloudMat(replaceLinInd-jj);
        else
            replaceLinInd1=replaceLinInd;
        end
    end

    [~,cInd]=ind2sub(size(data.DBZ),replaceLinInd1);
    gasAttCloud(cInd)=gasAttCloud(cInd)'+gasAttCloudMat(replaceLinInd1);

    %PIA gas
    piaGasMat2=cumsum(gasAttCloudMat,1).*2;

    data.DBZcorrGas=data.DBZ+piaGasMat2;

    %% Remove all up pointing and unsuitable data

    dbzOrig=data.DBZ;

    % Noise source cal (10), missing (11)
    badInds=find(any(data.FLAG>9,1));

    % Land surface
    badInds=cat(2,badInds,find(any(data.FLAG==8,1)));

    % Zenith (2), pointing (3), scanning (4), transision (5), failure (6)
    badInds=cat(2,badInds,find(data.ANTFLAG>1));

    infields=fields(data);
    for bb=1:length(infields)
        if strcmp(infields{bb},'DBZ') | strcmp(infields{bb},'FLAG') | ...
                strcmp(infields{bb},'rotation') | strcmp(infields{bb},'elevation')
            currfield=data.(infields{bb});
            currfield(:,badInds)=nan;
            data.(infields{bb})=currfield;
        end
    end

    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG~=1)=nan;

    data.dbzMaskedCorrGas=data.DBZcorrGas;
    data.dbzMaskedCorrGas(data.FLAG~=1)=nan;

    %% Get hydrometeor attenuation

    [sig0measAtt,surfFlag,refSig0,sig0model,piaHydromet2]=getRefAtten_fromGasCorr(data);

    %% Hitschfeld Bordan from surface up
    piaHydromet1=piaHydromet2/2;

    % Fix PIA
    piaHydrometInt=piaHydromet1;
    piaHydrometInt(piaHydrometInt<0)=0;

    % Handle extinct
    extInds=any(data.FLAG==3,1);
    extInds=imdilate(extInds,strel('disk',3));
    piaHydrometInt(isnan(piaHydrometInt) & extInds==0)=0;

    piaHydrometInt=fillmissing(piaHydrometInt,'linear','EndValues','nearest');

    cloudInds=any(data.FLAG==1,1);
    piaHydrometInt(~cloudInds)=nan;
    piaHydrometInt(badInds)=nan;

    zHB=hitschfeldBordan_surfUp(data.dbzMaskedCorrGas,piaHydrometInt,data.range);

    % Prepare plot
    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);

    refSig0(badInds)=nan;

    %% Plot in 30 min increments

    disp('Plotting ...');

    startPlot=startTime;

    while startPlot<endTime

        close all

        endPlot=startPlot+minutes(30);
        timeInds=find(data.time>=startPlot & data.time<=endPlot);

        if sum(sum(~isnan(zHB(:,timeInds))))==0
            startPlot=endPlot;
            continue
        end

        close all

        f1 = figure('Position',[200 500 1800 1250],'DefaultAxesFontSize',12,'renderer','painters','visible',showPlot);

        colormap jet

        s1=subplot(4,1,1);
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,dbzOrig(:,timeInds),'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([-0.1 upperLimDBZ]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Reflectivity (dBZ)')

        s2=subplot(4,1,2);
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.dbzMaskedCorrGas(:,timeInds),'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([-0.1 upperLimDBZ]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Reflectivity corrected for gaseous attenuation (dBZ)')

        s3=subplot(4,1,3);
        hold on
        surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,zHB(:,timeInds),'edgecolor','none');
        view(2);
        ylabel('Altitude (km)');
        caxis([-25 25]);
        ylim([-0.1 upperLimDBZ]);
        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
        colorbar
        grid on
        title('Reflectivity corrected for gaseous and liquid attenuation (dBZ)')

        s4=subplot(4,1,4);
        hold on
        l0=plot(data.time(timeInds),sig0model(timeInds),'-c','linewidth',2);
        l1=plot(data.time(timeInds),sig0measClear(timeInds),'-b','linewidth',1);
        l2=plot(data.time(timeInds),sig0measCloud(timeInds),'color',[0.5 0.5 0.5],'linewidth',0.5);
        l3=plot(data.time(timeInds),refSig0(timeInds),'-r','linewidth',2);
        ylabel('Sig0 (dB)');
        ylim([0 20]);

        yyaxis right
        l6=plot(data.time(timeInds),gasAttCloud(timeInds)*2,'-k','linewidth',1);
        l7=plot(data.time(timeInds),piaHydrometInt(timeInds)*2,'-','color',[0,0.5,0],'linewidth',2);
        l8=plot(data.time(timeInds),piaHydromet2(timeInds),'-g','linewidth',1);
        ylabel('Atten. (dB)');
        ylim([-5 15]);
        grid on
        set(gca,'YColor','k');

        xlim([data.time(timeInds(1)),data.time(timeInds(end))]);

        legend([l0 l1 l2 l3 l6 l7 l8],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
            'sig0 ref','2-way gas att','2-way liq att int','2-way liq att'},...
            'orientation','horizontal','location','north');
        title([datestr(data.time(timeInds(1)),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(timeInds(end)),'yyyy-mm-dd HH:MM:SS')])

        s1Pos=s1.Position;
        s2Pos=s2.Position;
        s3Pos=s3.Position;
        s4Pos=s4.Position;

        s1.Position=[s1Pos(1),s1Pos(2),s4Pos(3),s1Pos(4)];
        s2.Position=[s2Pos(1),s2Pos(2),s4Pos(3),s2Pos(4)];
        s3.Position=[s3Pos(1),s3Pos(2),s4Pos(3),s3Pos(4)];

        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_Flight',num2str(aa),'_att_',datestr(data.time(timeInds(1)),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(timeInds(end)),'yyyymmdd_HHMMSS')],'-dpng','-r0')

        startPlot=endPlot;
    end
end