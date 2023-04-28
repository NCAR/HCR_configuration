% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='cset'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';

showPlot='on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'attenPlots/cases/'];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/atten_',project,'.txt'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)

    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);

    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    %% Get data

    disp("Getting data ...");

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
    [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_gas_atten(data);
    piaGasMat2=cumsum(gasAttCloudMat,1).*2;

    data.DBZcorrGas=data.DBZ+piaGasMat2;

    %% Remove all up pointing and unsuitable data

    dbzOrig=data.DBZ;

    % Noise source cal (10), missing (11)
    badInds=find(any(data.FLAG>9,1));
    % Extinct
    % badInds=cat(2,badInds,find(any(data.FLAG==3,1)));
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

   %% Get reference attenuation
   [sig0measAtt,surfFlag,refSig0,refFlag,sig0model,piaHydromet2]=getRefAtten_fromGasCorr(data);
   
   %% Hitschfeld Bordan from surface up
   piaHydromet1=piaHydromet2/2;
   zHB=hitschfeldBordan_surfUp(data.dbzMaskedCorrGas,piaHydromet1,data.range);

    %% Plot lines

    close all

    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);

    refSig0(badInds)=nan;
    
    sig0refMeas=nan(size(data.time));
    sig0refMeas(refFlag==1)=refSig0(refFlag==1);
    sig0refInt=nan(size(data.time));
    sig0refInt(refFlag==2)=refSig0(refFlag==2);
    sig0refMod=nan(size(data.time));
    sig0refMod(refFlag==3)=refSig0(refFlag==3);

    f1 = figure('Position',[200 500 1800 500],'DefaultAxesFontSize',12,'renderer','painters');

    hold on
    l0=plot(data.time,sig0model,'-c','linewidth',2);
    l1=plot(data.time,sig0measClear,'-b','linewidth',1);
    l2=plot(data.time,sig0measCloud,'color',[0.5 0.5 0.5],'linewidth',0.5);
    l3=plot(data.time,sig0refMeas,'-r','linewidth',2);
    l4=plot(data.time,sig0refInt,'-','color',[0.5 0 1],'linewidth',2);
    l5=plot(data.time,sig0refMod,'-m','linewidth',2);
    ylabel('Sig0 (dB)');
    ylim([0 20]);

    yyaxis right
    l6=plot(data.time,gasAttCloud*2,'-k','linewidth',1);
    l7=plot(data.time,piaHydromet2,'-g','linewidth',1);
    ylabel('Atten. (dB)');
    ylim([-5 15]);
    grid on
    set(gca,'YColor','k');

    xlim([data.time(1),data.time(end)]);

    legend([l0 l1 l2 l3 l4 l5 l6 l7],{'sig0 model','sig0 meas clear','sig0 meas cloud',...
        'sig0 ref meas','sig0 ref int','sig0 ref mod','2-way gas att','2-way liq att'},...
        'orientation','horizontal','location','north');
    title([datestr(data.time(1),'yyyy-mm-dd HH:MM:SS'),' to ',datestr(data.time(end),'yyyy-mm-dd HH:MM:SS')])

    set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_sig0att_',...
            datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
    %% Plot refl

    timeInds=1:round(length(data.time)/3000):length(data.time);

    ylimRefl=ceil(max(data.asl(~isnan(zHB)))./1000);

    f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'renderer','painters');

    colormap jet

    subplot(3,1,1)
    hold on
    surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,dbzOrig(:,timeInds),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-25 25]);
    ylim([-0.1 ylimRefl]);
    xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
    colorbar
    grid on
    title('Reflectivity (dBZ)')

    subplot(3,1,2)
    hold on
    surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,data.dbzMaskedCorrGas(:,timeInds),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-25 25]);
    ylim([-0.1 ylimRefl]);
    xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
    colorbar
    grid on
    title('Reflectivity corrected for gaseous attenuation (dBZ)')

    subplot(3,1,3)
    hold on
    surf(data.time(:,timeInds),data.asl(:,timeInds)./1000,zHB(:,timeInds),'edgecolor','none');
    view(2);
    ylabel('Altitude (km)');
    caxis([-25 25]);
    ylim([-0.1 ylimRefl]);
    xlim([data.time(timeInds(1)),data.time(timeInds(end))]);
    colorbar
    grid on
    title('Reflectivity corrected for gaseous and liquid attenuation (dBZ)')

    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_dbz_',...
        datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
end