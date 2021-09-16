% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

figdir=['/scr/sleet2/rsfdata/projects/spicule/hcr/',quality,'/cfradial/',qcVersion,'/deAliasVEL/'];

ylimits=[0 13];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velDeAlias_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=3:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    
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
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    if isfield(data,'VEL_CORR')
        data.VEL_MASKED=data.VEL_CORR;
        data=rmfield(data,'VEL_CORR');
    end
    
    %% Correct velocity folding
    
    disp('De-aliasing ...');
    velDeAliased1=dealiasAreaPos(data.VEL_MASKED,data.elevation);   
    velDeAliased=dealiasAreaNeg(velDeAliased1,data.elevation);   
    
    %% Plot
    
    disp('Plotting ...');
    
    close all
    
    f1=figure('DefaultAxesFontSize',12,'Position',[0 300 1700 1200],'visible','on');
    
    s1=subplot(3,1,1);
    surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-8 8]);
    colormap(s1,jet);
    colorbar;
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
    
    s2=subplot(3,1,2);
    surf(data.time,data.asl./1000,velDeAliased1,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-15 15]);
    colormap(s2,jet);
    colorbar;
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
    
    s3=subplot(3,1,3);
    surf(data.time,data.asl./1000,velDeAliased,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-15 15]);
    colormap(s3,jet);
    colorbar;
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
    
    linkaxes([s1 s2 s3],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velDeAlias_',...
        datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
    
    
    
end