% Analyze HCR clouds

clear all;
close all;

plotTest=1;

project='otrec'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz

if strcmp(project,'otrec')
    ylimits=[-0.2 15];
else
    ylimits=[-0.2 10];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

%figdir=['/h/eol/romatsch/hcrCalib/clouds/cloudPuzzle/'];
figdir=['/home/romatsch/plots/HCR/cloudPuzzle/',project,'/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/10hz/'];

% Loop through cases
casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/cloudPuzzle_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
    data.DBZ=[];
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
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
    
     %% Cloud puzzle
    disp('Making cloud puzzle');
    
    cloudPuzzle=f_cloudPuzzle_radial(data);
    uClouds=unique(cloudPuzzle(~isnan(cloudPuzzle)));
    cloudCount=length(uClouds);
    

    %% Plot
    
    disp('Plotting ...');
    
    close all
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,100,1300,900]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% DBZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1=subplot(2,1,1);
    hold on;
    
    sub1=surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity')
    grid on
    
    ax2=subplot(2,1,2);
    
    colMapIn=jet(cloudCount-1);
    % Make order random
    indsCol=randperm(size(colMapIn,1));
    colMapInds=cat(2,indsCol',colMapIn);
    colMapInds=sortrows(colMapInds);
    colMap=cat(1,[0 0 0],colMapInds(:,2:end));
    
    hold on;
    sub2=surf(data.time,data.asl./1000,cloudPuzzle,'edgecolor','none');
    view(2);
    ax2.Colormap=colMap;
    caxis([-0.5 cloudCount-1+0.5])
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Cloud Puzzle')
    grid on
    colorbar
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,project,'_',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_cloudPuzzle.png'],'-dpng','-r0');
    
end