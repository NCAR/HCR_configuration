% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
qcVersion='v1.0';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

ylimits=[0 13];

casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/velDeAlias_',project,'.txt'];

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'deAliasVEL/cases/'];

% Loop through cases

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=5:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    %% Load data
    
    disp('Loading data ...');
    
    data=[];
    data.nyquist_velocity=[];
    data.VEL_CORR=[];
    data.FLAG=[];
    data.ANTFLAG=[];
               
    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if length(fileList)==0
        disp('No data files found.');
        return
    end
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    velMasked=data.VEL_CORR;
    velMasked(:,data.ANTFLAG>2)=nan;
    velMasked(data.FLAG~=1)=nan;
      
    %% Correct velocity folding
    
    disp('De-aliasing ...');
    velDeAliased=dealiasArea(velMasked,data.elevation,data.nyquist_velocity);
       
    %% Plot
    
    disp('Plotting ...');
    
    close all
    
    f1=figure('DefaultAxesFontSize',12,'Position',[0 300 1700 1200],'visible','on');
    
    s1=subplot(2,1,1);
    surf(data.time,data.asl./1000,velMasked,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-8 8]);
    colormap(s1,jet);
    colorbar;
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
    
    s2=subplot(2,1,2);
    surf(data.time,data.asl./1000,velDeAliased,'edgecolor','none');
    view(2);
    ylim(ylimits);
    xlim([data.time(1),data.time(end)]);
    caxis([-15 15]);
    colormap(s2,jet);
    colorbar;
    ylabel('Altitude (km)');
    title(['HCR radial velocity (m s^{-1})']);
    grid on
        
    linkaxes([s1 s2],'xy');
    
    set(gcf,'PaperPositionMode','auto')
    print(f1,[figdir,project,'_velDeAlias_',...
        datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
       
end