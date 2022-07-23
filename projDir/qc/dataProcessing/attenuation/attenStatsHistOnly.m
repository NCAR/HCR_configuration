% Calculate liquid water content from HCR ocean return

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='otrec'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.1';

showPlot='on';
ylimRefl=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

dataDir=HCRdir(project,quality,qcVersion,freqData);

figdir=[dataDir(1:end-5),'attenPlots/stats/cases/'];

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

% Inizialize
dbzCY=[];
dbzCN=[];
dbzWY=[];
dbzWN=[];

ldrCY=[];
ldrCN=[];
ldrWY=[];
ldrWN=[];

velCY=[];
velCN=[];
velWY=[];
velWN=[];

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    %% Get data

    disp("Getting data ...");

    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];

    data.DBZ = [];
    data.LDR=[];
    data.VEL_MASKED=[];
    data.U_SURF=[];
    data.V_SURF=[];
    data.SST=[];
    data.TEMP=[];
    data.PRESS=[];
    data.RH=[];
    data.TOPO=[];
    data.FLAG=[];
    data.rotation=[];
    data.MELTING_LAYER=[];
    data.pulse_width=[];

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

    data.frq=ncread(fileList{1},'frequency');

    %% Remove all up pointing data

    dbzOrig=data.DBZ;

    upInds=find(data.elevation>-85);
    upInds=cat(2,upInds,find(any(data.FLAG>9,1)),find(any(data.FLAG==3,1)),find(any(data.FLAG==8,1)));

    infields=fields(data);
    for bb=1:length(infields)
        if strcmp(infields{bb},'DBZ') | strcmp(infields{bb},'FLAG') | ...
                strcmp(infields{bb},'rotation') | strcmp(infields{bb},'elevation')
            currfield=data.(infields{bb});
            currfield(:,upInds)=nan;
            data.(infields{bb})=currfield;
        end
    end

    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;

    data.ldrMasked=data.LDR;
    data.ldrMasked(data.FLAG>1)=nan;

    meltMask=data.MELTING_LAYER;
    meltMask(data.FLAG~=1)=nan;

   %% Get reference attenuation
   [sig0measAtt,surfFlag,refSig0,refFlag,sig0model,piaGas2,piaHydromet2]=getRefAtten(data);

   %% Divede into categories

   % Flag: 1 attenuation, 0 no attenuation, nan no cloud or extinct
   ynMat=nan(size(data.DBZ));
   ynMat(:,piaHydromet2>=1)=1;
   ynMat(:,piaHydromet2<1)=0;

   % Find where only warm or only cold clouds   
   cYes=any(meltMask>=20,1);
   wYes=any(meltMask<20,1);

   cYesMat=repmat(cYes,size(data.DBZ,1),1);
   wYesMat=repmat(wYes,size(data.DBZ,1),1);

   % Indices for cold yes, cold no, warm yes, warm no
   cy=find(ynMat==1 & cYesMat==1 & wYesMat==0);
   cn=find(ynMat==0 & cYesMat==1 & wYesMat==0);
   wy=find(ynMat==1 & cYesMat==0 & wYesMat==1);
   wn=find(ynMat==0 & cYesMat==0 & wYesMat==1);

   % DBZ
   dbzCY=cat(1,dbzCY,data.dbzMasked(cy));
   dbzCN=cat(1,dbzCN,data.dbzMasked(cn));
   dbzWY=cat(1,dbzWY,data.dbzMasked(wy));
   dbzWN=cat(1,dbzWN,data.dbzMasked(wn));

   % LDR
   ldrCY=cat(1,ldrCY,data.ldrMasked(cy));
   ldrCN=cat(1,ldrCN,data.ldrMasked(cn));
   ldrWY=cat(1,ldrWY,data.ldrMasked(wy));
   ldrWN=cat(1,ldrWN,data.ldrMasked(wn));

   % VEL
   velCY=cat(1,velCY,data.VEL_MASKED(cy));
   velCN=cat(1,velCN,data.VEL_MASKED(cn));
   velWY=cat(1,velWY,data.VEL_MASKED(wy));
   velWN=cat(1,velWN,data.VEL_MASKED(wn));

    %% Plot lines

    close all

    sig0measClear=nan(size(data.time));
    sig0measClear(surfFlag==2)=sig0measAtt(surfFlag==2);
    sig0measCloud=nan(size(data.time));
    sig0measCloud(surfFlag==1)=sig0measAtt(surfFlag==1);

    refSig0(upInds)=nan;
    piaGas2(upInds)=nan;

    sig0refMeas=nan(size(data.time));
    sig0refMeas(refFlag==1)=refSig0(refFlag==1);
    sig0refInt=nan(size(data.time));
    sig0refInt(refFlag==2)=refSig0(refFlag==2);
    sig0refMod=nan(size(data.time));
    sig0refMod(refFlag==3)=refSig0(refFlag==3);

    f1 = figure('Position',[200 700 1800 500],'DefaultAxesFontSize',12,'renderer','painters');

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
    l6=plot(data.time,piaGas2,'-k','linewidth',1);
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
                
end

%% Histograms counts

edgesDBZ=-60:5:30;
edgesLDR=-30:2:10;
edgesVEL=-10:10;
    
numDBZcy=histcounts(dbzCY,edgesDBZ);
numDBZcn=histcounts(dbzCN,edgesDBZ);
numDBZwy=histcounts(dbzWY,edgesDBZ);
numDBZwn=histcounts(dbzWN,edgesDBZ);

numLDRcy=histcounts(ldrCY,edgesLDR);
numLDRcn=histcounts(ldrCN,edgesLDR);
numLDRwy=histcounts(ldrWY,edgesLDR);
numLDRwn=histcounts(ldrWN,edgesLDR);

numVELcy=histcounts(velCY,edgesVEL);
numVELcn=histcounts(velCN,edgesVEL);
numVELwy=histcounts(velWY,edgesVEL);
numVELwn=histcounts(velWN,edgesVEL);

%% Plot histograms

close all

f1 = figure('Position',[200 700 900 1000],'DefaultAxesFontSize',12);
colmapC=[1,0,1;0,1,1];
colmapW=[1,0,0;0,0,1];

% DBZ
subplot(3,2,1);

b=bar(edgesDBZ(1:end-1)+(edgesDBZ(2)-edgesDBZ(1))/2,(cat(1,numDBZcy./sum(numDBZcy),numDBZcn./sum(numDBZcn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapC(k,:);
end
legend('Att','No att','Location','northwest')
xlabel('DBZ (dBZ)');
ylabel('Percent (%)');
title('Cold');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

subplot(3,2,2);

b=bar(edgesDBZ(1:end-1)+(edgesDBZ(2)-edgesDBZ(1))/2,(cat(1,numDBZwy./sum(numDBZwy),numDBZwn./sum(numDBZwn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapW(k,:);
end
legend('Att','No att','Location','northwest')
xlabel('DBZ (dBZ)');
ylabel('Percent (%)');
title('Warm');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

% LDR

subplot(3,2,3);

b=bar(edgesLDR(1:end-1)+(edgesLDR(2)-edgesLDR(1))/2,(cat(1,numLDRcy./sum(numLDRcy),numLDRcn./sum(numLDRcn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapC(k,:);
end
xlabel('LDR (dB)');
ylabel('Percent (%)');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

subplot(3,2,4);

b=bar(edgesLDR(1:end-1)+(edgesLDR(2)-edgesLDR(1))/2,(cat(1,numLDRwy./sum(numLDRwy),numLDRwn./sum(numLDRwn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapW(k,:);
end
xlabel('LDR (dB)');
ylabel('Percent (%)');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

% VEL

subplot(3,2,5);

b=bar(edgesVEL(1:end-1)+(edgesVEL(2)-edgesVEL(1))/2,(cat(1,numVELcy./sum(numVELcy),numVELcn./sum(numVELcn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapC(k,:);
end
xlabel('VEL (dB)');
ylabel('Percent (%)');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

subplot(3,2,6);

b=bar(edgesVEL(1:end-1)+(edgesVEL(2)-edgesVEL(1))/2,(cat(1,numVELwy./sum(numVELwy),numVELwn./sum(numVELwn)))'.*100,1,'FaceColor','flat');
for k=1:2
    b(k).CData=colmapW(k,:);
end
xlabel('VEL (dB)');
ylabel('Percent (%)');
set(gca, 'XGrid', 'on', 'YGrid', 'off')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_histogram'],'-dpng','-r0')