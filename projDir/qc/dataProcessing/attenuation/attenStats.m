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
dbzC=[];
dbzW=[];
ldrC=[];
ldrW=[];
velC=[];
velW=[];
attC=[];
attW=[];

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

   piaHydromet2(piaHydromet2>20)=nan;
   piaHydromet2(piaHydromet2<0)=nan;
   %% Divede into only cold and only warm categories

   % Find where only warm or only cold clouds   
   cYes=any(meltMask>=20,1);
   wYes=any(meltMask<20,1);

   cInds=find(cYes==1 & wYes==0 & ~isnan(piaHydromet2));
   wInds=find(cYes==0 & wYes==1 & ~isnan(piaHydromet2));

   dbzC=cat(2,dbzC,data.dbzMasked(:,cInds));
   dbzW=cat(2,dbzW,data.dbzMasked(:,wInds));
   ldrC=cat(2,ldrC,data.ldrMasked(:,cInds));
   ldrW=cat(2,ldrW,data.ldrMasked(:,wInds));
   velC=cat(2,velC,data.VEL_MASKED(:,cInds));
   velW=cat(2,velW,data.VEL_MASKED(:,wInds));
   attC=cat(2,attC,piaHydromet2(cInds));
   attW=cat(2,attW,piaHydromet2(wInds));
   
end

%% Cloud pixel sum

pixSumC=sum(~isnan(dbzC),1);
pixSumW=sum(~isnan(dbzW),1);

close all

f1 = figure('Position',[200 700 600 1000],'DefaultAxesFontSize',12);

subplot(2,1,1)
scatter(pixSumC,attC,'filled');
title('Cold cloud pixels vs hydrometeor attenuation')
xlabel('Cloud pixel number')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,1,2)
scatter(pixSumW,attW,'filled');
title('Warm cloud pixels vs hydrometeor attenuation')
xlabel('Cloud pixel number')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_pixNumVSatt'],'-dpng','-r0')

%% Max refl and refl sum

maxReflC=max(dbzC,[],1,'omitnan');
maxReflW=max(dbzW,[],1,'omitnan');

linDBZc=10.^(dbzC./10);
linDBZw=10.^(dbzW./10);

reflSumClin=sum(linDBZc,1,'omitnan');
reflSumWlin=sum(linDBZw,1,'omitnan');

close all

f1 = figure('Position',[200 700 1300 1000],'DefaultAxesFontSize',12);

subplot(2,2,1)
scatter(maxReflC,attC,'filled');
title('Cold max refl. vs hydrometeor attenuation')
xlabel('Max reflectivity (dBZ)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,2)
scatter(maxReflW,attW,'filled');
title('Warm max refl. vs hydrometeor attenuation')
xlabel('Max reflectivity (dBZ)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,3)
scatter(reflSumClin,attC,'filled');
title('Cold refl. sum vs hydrometeor attenuation')
xlabel('Reflectivity sum (Z)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,4)
scatter(reflSumWlin,attW,'filled');
title('Warm refl. sum vs hydrometeor attenuation')
xlabel('Reflectivity sum (Z)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_reflVSatt'],'-dpng','-r0')

