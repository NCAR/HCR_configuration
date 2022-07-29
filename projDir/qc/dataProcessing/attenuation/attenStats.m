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
    data.ICING_LEVEL=[];
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

    %% Get reference attenuation
    [sig0measAtt,surfFlag,refSig0,refFlag,sig0model,piaGas2,piaHydromet2]=getRefAtten(data);

    piaHydromet2(piaHydromet2>20)=nan;
    piaHydromet2(piaHydromet2<0)=nan;

    %% Divede into only cold and only warm categories

    % Find where only warm or only cold clouds
    meltMask=data.MELTING_LAYER;
    meltMask(data.asl<data.ICING_LEVEL+200 & data.asl>data.ICING_LEVEL-200)=50;
    meltMask(data.FLAG~=1)=nan;

    cYes=any(meltMask>=20,1);
    mYes=any(meltMask==50,1);
    wYes=any(meltMask<20,1);

    cInds=find(cYes==1 & wYes==0 & mYes==0 & ~isnan(piaHydromet2));
    wInds=find(cYes==0 & wYes==1 & mYes==0 & ~isnan(piaHydromet2));

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
colormap('jet');

subplot(2,1,1)
centers={5:10:ceil(max(pixSumC)/10)*10-5 0.25:0.5:ceil(max(attC))-0.25};

xTickLoc=0.5:5:length(centers{1})+0.5;
xTickLabel=num2str(((0:5:length(centers{1}))*10)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,pixSumC',attC');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 length(centers{1})+0.5])
ylim([0.5 length(centers{2})+0.5])

title('Cold cloud pixels vs hydrometeor attenuation')
xlabel('Cloud pixel number')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,1,2)
centers={5:10:ceil(max(pixSumW)/10)*10-5 0.25:0.5:ceil(max(attW))-0.25};

xTickLoc=0.5:5:length(centers{1})+0.5;
xTickLabel=num2str(((0:5:length(centers{1}))*10)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,pixSumW',attW');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 length(centers{1})+0.5])
ylim([0.5 length(centers{2})+0.5])

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

reflSumC=10.*log10(reflSumClin);
reflSumW=10.*log10(reflSumWlin);

f1 = figure('Position',[200 700 1300 1000],'DefaultAxesFontSize',12);
colormap('jet');

subplot(2,2,1)
centers={-57.5:2.5:25 0.25:0.5:ceil(max(attC))-0.25};

xTickLoc=0.5:4:50;
xTickLabel=num2str((-60:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,maxReflC',attC');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 35.5])
ylim([0.5 length(centers{2})+0.5])

title('Cold max refl. vs hydrometeor attenuation')
xlabel('Max reflectivity (dBZ)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,2)
centers={-57.5:2.5:25 0.25:0.5:ceil(max(attW))-0.25};

xTickLoc=0.5:4:50;
xTickLabel=num2str((-60:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,maxReflW',attW');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 35.5])
ylim([0.5 length(centers{2})+0.5])
title('Warm max refl. vs hydrometeor attenuation')
xlabel('Max reflectivity (dBZ)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,3)
centers={-57.5:2.5:40 0.25:0.5:ceil(max(attC))-0.25};

xTickLoc=0.5:4:50;
xTickLabel=num2str((-60:10:40)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,reflSumC',attC');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 40.5])
ylim([0.5 length(centers{2})+0.5])

title('Cold refl. sum vs hydrometeor attenuation')
xlabel('Reflectivity sum (Z)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,4)

centers={-57.5:2.5:40 0.25:0.5:ceil(max(attW))-0.25};

xTickLoc=0.5:4:50;
xTickLabel=num2str((-60:10:40)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,reflSumW',attW');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 40.5])
ylim([0.5 length(centers{2})+0.5])

title('Warm refl. sum vs hydrometeor attenuation')
xlabel('Reflectivity sum (Z)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_reflVSatt'],'-dpng','-r0')

%% Min ldr and med ldr

minLdrC=min(ldrC,[],1,'omitnan');
minLdrW=min(ldrW,[],1,'omitnan');

ldrMedC=median(ldrC,1,'omitnan');
ldrMedW=median(ldrW,1,'omitnan');

f1 = figure('Position',[200 700 1300 1000],'DefaultAxesFontSize',12);
colormap('jet');

subplot(2,2,1)
centers={-32.5:2.5:25 0.25:0.5:ceil(max(attC))-0.25};

xTickLoc=2.5:4:50;
xTickLabel=num2str((-30:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,minLdrC',attC');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 24.5])
ylim([0.5 length(centers{2})+0.5])

title('Cold min LDR vs hydrometeor attenuation')
xlabel('Min LDR (dB)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,2)

centers={-32.5:2.5:25 0.25:0.5:ceil(max(attW))-0.25};

xTickLoc=2.5:4:50;
xTickLabel=num2str((-30:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,minLdrW',attW');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 16.5])
ylim([0.5 length(centers{2})+0.5])

title('Warm min LDR vs hydrometeor attenuation')
xlabel('Min LDR (dB)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,3)
centers={-32.5:2.5:25 0.25:0.5:ceil(max(attC))-0.25};

xTickLoc=2.5:4:50;
xTickLabel=num2str((-30:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,ldrMedC',attC');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 24.5])
ylim([0.5 length(centers{2})+0.5])

title('Cold LDR median vs hydrometeor attenuation')
xlabel('LDR (dB)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,2,4)

centers={-32.5:2.5:25 0.25:0.5:ceil(max(attW))-0.25};

xTickLoc=2.5:4:50;
xTickLabel=num2str((-30:10:30)');
yTickLoc=0.5:4:length(centers{2})+0.5;
yTickLabel=num2str((0:2:50)');

xy=cat(2,ldrMedW',attW');
xy(any(isnan(xy),2),:)=[];

N=hist3(xy,'Ctrs',centers);
N(N==0)=nan;

imagesc(N','AlphaData',~isnan(N'))
set(gca,'YDir','normal');
set(gca,'Xtick',xTickLoc);
set(gca,'Xticklabel',xTickLabel);
set(gca,'Ytick',yTickLoc);
set(gca,'Yticklabel',yTickLabel);

caxis([0 max(max(N))/4])
colorbar

xlim([0.5 16.5])
ylim([0.5 length(centers{2})+0.5])

title('Warm LDR median vs hydrometeor attenuation')
xlabel('LDR (dB)')
ylabel('Hydrometeor 2 way attenuation (dB)')
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_ldrVSatt'],'-dpng','-r0')

%% No LDR

noLdrC=find(sum(~isnan(ldrC),1)==0);
noLdrW=find(sum(~isnan(ldrW),1)==0);

attNoLdrC=attC(noLdrC);
attNoLdrW=attW(noLdrW);

f1 = figure('Position',[200 700 600 1000],'DefaultAxesFontSize',12);

subplot(2,1,1)
edgesNoLdr=0:0.1:10;

hC=histcounts(attNoLdrC,edgesNoLdr);

bar(edgesNoLdr(1:end-1)+0.05,hC,1)

xlim([0 10])

title('Hydrometeor attenuation for rays with no LDR, cold')
ylabel('Ray counts')
xlabel('Hydrometeor 2 way attenuation (dB)')
grid on

subplot(2,1,2)

hW=histcounts(attNoLdrW,edgesNoLdr);

bar(edgesNoLdr(1:end-1)+0.05,hW,1)

xlim([0 10])

title('Hydrometeor attenuation for rays with no LDR, warm')
ylabel('Ray counts')
xlabel('Hydrometeor 2 way attenuation (dB)')
grid on

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_noLDRatt'],'-dpng','-r0')