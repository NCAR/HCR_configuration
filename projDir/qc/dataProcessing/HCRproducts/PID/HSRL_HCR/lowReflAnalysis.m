% Plot HCR pid from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='combined'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

showPlot='off';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-4),'pidPlotsComb/lowRefl/'];

% cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
% units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
%     'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

dbzAll=[];
pidAll=[];
dbzAllNP=[];
pidAllNP=[];

for aa=1:15
    disp(['Flight ',num2str(aa)]);
    disp('Loading data ...')


    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    %% Get HCR data

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    data=[];
    data.HCR_DBZ=[];
    data.PID=[];
   
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

    if isempty(data.time)
        startTime=endTime;
        continue
    end

    DBZorig=data.HCR_DBZ;

    %% All

    data.HCR_DBZ(data.HCR_DBZ>-10)=nan;
    data.HCR_DBZ(data.PID>9)=nan;
    data.PID(isnan(data.HCR_DBZ))=nan;

    dbzAll=cat(1,dbzAll,data.HCR_DBZ(~isnan(data.HCR_DBZ)));
    pidAll=cat(1,pidAll,data.PID(~isnan(data.HCR_DBZ)));

    %% Near plane gates

    DBZNP=data.HCR_DBZ;

    for kk=1:size(DBZorig,2)
        dbzCol=DBZorig(:,kk);
        firstHigh=min(find(dbzCol>0));
        if ~isempty(firstHigh)
            DBZNP(firstHigh:end,kk)=nan;
        end
    end

    data.PID(isnan(DBZNP))=nan;

    dbzAllNP=cat(1,dbzAllNP,DBZNP(~isnan(DBZNP)));
    pidAllNP=cat(1,pidAllNP,data.PID(~isnan(DBZNP)));

end

%% Plot histogram
close all

units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
    'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};

f1=figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible','on');

scatter(pidAll,dbzAll);

ylabel('Reflectivity (dBZ)');

xlim([0.5,9.5]);
xticklabels(units_str_hcr);

title([project,' categories in low reflectivity regions']);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_lowDBZ.png'],'-dpng','-r0')

save([figdir,'pid_dbzLow.mat'],'dbzAll','pidAll');

%% Below -30

close all

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];

dbzPID=cat(2,dbzAll,pidAll);
dbzPID(find(dbzAll>-30),:)=[];

edges=0.5:9.5;
countsCat=histcounts(dbzPID,edges);

f1=figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible','on');

b=bar(1:9,countsCat./sum(countsCat).*100,1);

b.FaceColor='flat';
for jj=1:9
    b.CData(jj,:)=cscale_hcr(jj,:);
end

xlim([0.5,9.5]);
ylim([0 50]);

xticklabels(units_str_hcr);

title([project,' percent per category DBZ<=-30 dBZ']);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

sumLiq=sum(countsCat(5:6));
sumFrozen=countsCat(9);

cloudPerc=sumLiq./sum(countsCat).*100;
sFr=sumFrozen./sum(countsCat).*100;

text(1,46,['Cloud liquid: ',num2str(cloudPerc,2),'%. Small frozen: ',num2str(sFr,2),'%.'],'fontsize',12,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_catsBelowMinus30dBZ.png'],'-dpng','-r0')

%% Below -30

close all

cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];

dbzPIDNP=cat(2,dbzAllNP,pidAllNP);
dbzPIDNP(find(dbzAllNP>-30),:)=[];

edges=0.5:9.5;
countsCat=histcounts(dbzPIDNP,edges);

f1=figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible','on');

b=bar(1:9,countsCat./sum(countsCat).*100,1);

b.FaceColor='flat';
for jj=1:9
    b.CData(jj,:)=cscale_hcr(jj,:);
end

xlim([0.5,9.5]);
ylim([0 50]);

xticklabels(units_str_hcr);

title([project,' percent per category DBZ<=-30 dBZ near aircraft pixels']);
set(gca, 'YGrid', 'on', 'XGrid', 'off')

sumLiq=sum(countsCat(5:6));
sumFrozen=countsCat(9);

cloudPerc=sumLiq./sum(countsCat).*100;
sFr=sumFrozen./sum(countsCat).*100;

text(1,46,['Cloud liquid: ',num2str(cloudPerc,2),'%. Small frozen: ',num2str(sFr,2),'%.'],'fontsize',12,'FontWeight','bold');

set(gcf,'PaperPositionMode','auto')
print(f1,[figdir,project,'_catsBelowMinus30dBZ_nearPlane.png'],'-dpng','-r0')

