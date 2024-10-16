% Analyze HCR clouds

clear all;
close all;

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, or 2hz

%unfoldVelocity=1;
offsetIn=-300;

thresholds.LDRlimits=[-16,-7]; % SOCRATES, OTREC, CSET default: [-16,-7]
thresholds.LDRspeclePix=50; % SOCRATES, OTREC, CSET default: [] (not used)
thresholds.LDRsolidity=0.5; % SOCRATES, OTREC, CSET default: [] (not used)
thresholds.LDRsearchPix=25; % SOCRATES, OTREC, CSET default: 18
thresholds.LDRstd=150; % SOCRATES, OTREC, CSET default: 100
thresholds.LDRaltDiff=100; % SOCRATES, OTREC, CSET default: 50

thresholds.VELsearchPix=80; % SOCRATES, OTREC, CSET default: 50
thresholds.VELstd=70; % SOCRATES, OTREC, CSET default: 35
thresholds.VELaltDiff=100; % SOCRATES, OTREC, CSET default: 100
thresholds.VELudDiff=-0.7; % SOCRATES, OTREC, CSET default: -0.7
thresholds.VEL_LDRdiff=600; % SOCRATES, OTREC, CSET default: 200

thresholds.outlier=350; % SOCRATES, OTREC, CSET default: 50
thresholds.length=10; % SOCRATES, OTREC, CSET default: 20

% Determines plot zoom.
if strcmp(project,'otrec')
    ylimits=[-0.2 15];
elseif strcmp(project,'socrates')
    ylimits=[-0.2 5];
elseif strcmp(project,'spicule')
    ylimits=[-0.2 13];
elseif strcmp(project,'cset')
    ylimits=[-0.2 9];
end

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/',quality,'/cfradial/',qcVersion,'/meltLayerPlots/cases/'];

if ~exist(figdir, 'dir')
    mkdir(figdir)
end

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/meltLayer_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);
%indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socrates/10hz/';

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=24:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
    data.DBZ=[];
    data.TEMP=[];
    data.WIDTH=[];
    data.FLAG=[];
    data.TOPO=[];
    
    dataVars=fieldnames(data);
    
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Check if VEL_MASKED is available
    try
        velTest=ncread(fileList{1},'VEL_MASKED');
        data.VEL_MASKED=[];
    catch
        data.VEL_CORR=[];
    end
    
    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
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
    
    if isfield(data,'VEL_MASKED')
        data.VEL_CORR=data.VEL_MASKED;
        data=rmfield(data,'VEL_MASKED');
    end
    
    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end
    
%     %% Correct velocity folding
%     
%     if unfoldVelocity
%         data.VEL_CORR=unfoldVel(data.VEL_CORR,data.FLAG,data.elevation);
%     end
    
    %% Find melting layer
    data.dbzMasked=data.DBZ;
    data.dbzMasked(data.FLAG>1)=nan;
    
    [meltLayer iceLayer offset]=f_meltLayer(data,offsetIn,thresholds);
    elevenInds=find(meltLayer==11);
    twelveInds=find(meltLayer==12);
    thirteenInds=find(meltLayer==13);
    fourteenInds=find(meltLayer==14);
    
    twentyoneInds=find(meltLayer==21);
    twentytwoInds=find(meltLayer==22);
    twentythreeInds=find(meltLayer==23);
    twentyfourInds=find(meltLayer==24);
    
    %% Plot
    
    timeMat=repmat(data.time,size(data.TEMP,1),1);
    ldrMasked=data.LDR;
    ldrMasked(data.FLAG>1)=nan;
    velMasked=data.VEL_CORR;
    velMasked(data.FLAG>1)=nan;
    
    close all
    
    if etime(datevec(endTime),datevec(startTime))<=900
        newInds=1:1:length(data.time);
    elseif etime(datevec(endTime),datevec(startTime))<=3600
        newInds=1:10:length(data.time);
    else
        newInds=1:100:length(data.time);
    end
    
    % Resample for plotting
    newDBZ=data.dbzMasked(:,newInds);
    newLDR=ldrMasked(:,newInds);
    newVEL=velMasked(:,newInds);
    newASL=data.asl(:,newInds);
    newTEMP=data.TEMP(:,newInds);
    newFindMelt=meltLayer(:,newInds);
    newTime=data.time(newInds);
    
    fig1=figure('DefaultAxesFontSize',11,'position',[100,1300,1500,1200],'visible','off');
    
    ax1=subplot(4,1,1);
    hold on;
    sub1=surf(newTime,newASL./1000,newDBZ,'edgecolor','none');
    view(2);
    sub1=colMapDBZ(sub1);
    scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
    scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
    scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
    scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
    
    scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
    scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');
    scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');
    scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Reflectivity and melting layer')
    grid on
    set(gca,'xticklabel',[])
    ax1.Position=[0.06 0.765 0.87 0.21];
    
    ax2=subplot(4,1,2);
    hold on;
    sub1=surf(newTime,newASL./1000,newFindMelt,'edgecolor','none');
    ax2.Colormap=([1 0 1;1 1 0]);
    view(2);
    scatter(timeMat(elevenInds),data.asl(elevenInds)./1000,10,'k','filled');
    scatter(timeMat(twelveInds),data.asl(twelveInds)./1000,10,'b','filled');
    scatter(timeMat(thirteenInds),data.asl(thirteenInds)./1000,10,'c','filled');
    scatter(timeMat(fourteenInds),data.asl(fourteenInds)./1000,10,'g','filled');
    
    scatter(timeMat(twentyoneInds),data.asl(twentyoneInds)./1000,10,'k','filled');
    scatter(timeMat(twentytwoInds),data.asl(twentytwoInds)./1000,10,'b','filled');
    scatter(timeMat(twentythreeInds),data.asl(twentythreeInds)./1000,10,'c','filled');
    scatter(timeMat(twentyfourInds),data.asl(twentyfourInds)./1000,10,'g','filled');
    
    plot(data.time,iceLayer./1000,'linewidth',1,'color',[0.6 0.6 0.6]);
    ax = gca;
    ax.SortMethod = 'childorder';
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('Melting layer')
    grid on
    set(gca,'xticklabel',[])
    ax2.Position=[0.06 0.525 0.87 0.21];
    
    %%%%%%%%%%%%%%%%%%%%%%%% LDR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax3=subplot(4,1,3);
    hold on;
    sub3=surf(newTime,newASL./1000,newLDR,'edgecolor','none');
    view(2);
    caxis([-25 5]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('LDR')
    grid on
    set(gca,'xticklabel',[])
    ax3.Position=[0.06 0.287 0.87 0.21];
    
    ax4=subplot(4,1,4);
    ax4.Colormap=jet;
    hold on;
    sub3=surf(newTime,newASL./1000,newVEL,'edgecolor','none');
    view(2);
    caxis([-8 8]);
    colorbar
    ylim(ylimits);
    ylabel('Altitude (km)');
    xlim([data.time(1),data.time(end)]);
    title('VEL')
    grid on
    ax4.Position=[0.06 0.05 0.87 0.21];
    
    linkaxes([ax1 ax2 ax3 ax4],'xy');
    
    formatOut = 'yyyymmdd_HHMM';
    set(gcf,'PaperPositionMode','auto')
    print([figdir,'meltLayer',datestr(startTime,formatOut),'_to_',datestr(endTime,formatOut)],'-dpng','-r0');
    
    if ~isempty(offset)
        disp(['Melting layer is on average ',num2str(offset),' m from the zero degree isotherm.'])
    end
end