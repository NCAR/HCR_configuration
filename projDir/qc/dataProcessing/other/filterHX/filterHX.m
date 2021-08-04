% Remove stripes in HX

clear all;
close all;

startTime=datetime(2021,6,21,3,50,0);
endTime=datetime(2021,6,21,3,56,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/v1.0/hxPlots/'];

% if ~exist(figdir, 'dir')
%     mkdir(figdir)
% end

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/other/filterHX/fixLDRtimes.txt'];

dataDir=HCRdir(project,quality,qcVersion,freqData);
%dataDir='/scr/sleet2/rsfdata/projects/spicule/hcr/qc1/cfradial/moments/10hz/';

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
    
    
    %% Load data
    data=[];
    data.DBMHX=[];
    data.LDR=[];
    data.FLAG=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(dataDir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
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
    
    %% Find indices where correction needs to be applied
    medianHXend=median(data.DBMHX(760:770,:));
    
    corrInd=medianHXend>-96;
    
    data.LDR(data.FLAG>1)=nan;
    
    dbm=data.DBMHX(:,corrInd==1);
    ldr=data.LDR(:,corrInd==1);
    
    %% Find filter indices
    med=medianHXend(corrInd==1);
    
    %dbmMask=dbm>-94;
    dbmMask=dbm>med+0.5;
           
    dbmMaskOut=zeros(size(dbmMask));
    
    for ii=1:size(dbmMask,2)
        rayIn=dbmMask(:,ii);
        dbmMaskOut(:,ii)=bwareaopen(rayIn,10);
    end
    
    maskUse=dbmMaskOut==0 & dbmMask==1;
    
    filterLDR=ldr;
    filterLDR(maskUse==1)=nan;
    
    LDR_masked=data.LDR;
    LDR_masked(:,corrInd==1)=filterLDR;
    
    maskUsePlot=zeros(size(data.LDR));
    maskUsePlot(:,corrInd==1)=maskUse;
    maskUsePlot(maskUsePlot==0)=nan;
    %% Plot
    ylimits=[-0.5 10];
    
    close all
    
    figure('DefaultAxesFontSize',11,'position',[1,100,1800,1200]);
    colormap jet
    
    ax1=subplot(3,1,1);
    fig1=surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
    view(2);
    ylim(ylimits);
    caxis([-40 -10])
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    
    ax1=subplot(3,1,2);
    fig1=surf(data.time,data.asl./1000,maskUsePlot,'edgecolor','none');
    view(2);
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Mask')
    
    ax1=subplot(3,1,3);
    fig1=surf(data.time,data.asl./1000,LDR_masked,'edgecolor','none');
    view(2);
    ylim(ylimits);
    caxis([-40 -10])
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR masked')
    %
    % linkaxes([ax1,ax2,ax3],'xy');
    
    formatOut = 'yyyymmdd_HHMM'; set(gcf,'PaperPositionMode','auto')
    print([figdir,datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut),'_ldrMask'],...
        '-dpng','-r0');
    
end