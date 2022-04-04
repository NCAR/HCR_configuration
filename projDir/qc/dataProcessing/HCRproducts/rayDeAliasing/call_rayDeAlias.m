% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
qcVersion='v1.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
dataFreq=10;

ylimits=[0 13];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velDeAlias_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'rayDeAlias/cases/'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    plotStart=datetime(2021,7,11,22,5,0);
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    data.nyquist_velocity=[];
    data.VEL_CORR=[];
    data.FLAG=[];
                   
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    velMasked=data.VEL_CORR;
    velMasked(data.FLAG~=1)=nan;
          
    %% Correct velocity folding
    
    nyq=mode(data.nyquist_velocity);
    
    disp('De-aliasing ...');
    velDeAliased=dealiasByRay(velMasked,data.elevation,nyq,dataFreq,data.time,plotStart);
           
    %% Plot
    
    disp('Plotting ...');
    
    close all
    
    f1=figure('DefaultAxesFontSize',12,'Position',[0 300 1700 1200],'visible','on');
    
    s1=subplot(2,1,1);
    h1=surf(data.time,data.asl./1000,velMasked,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-16 16]);
    colM=colormap(velCols);
    colormap(s1,colM);
    colorbar
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
    
    s2=subplot(2,1,2);
    h2=surf(data.time,data.asl./1000,velDeAliased,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-16 16]);
    colormap(s2,colM);
    colorbar
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
        
    linkaxes([s1 s2],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velDeAlias_',...
        datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
       
end