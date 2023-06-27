% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

plotIn.plotMR=0;
plotIn.plotMax=0;

convThresh=4;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-5),'pidPlots/cases/'];

if ~exist(figdir,'dir')
    mkdir(figdir)
end

% Loop through cases
casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/pid_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=1:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        %% Load data
        
        disp('Loading data');

        data=[];

        %HCR data
        data.DBZ_MASKED=[];
        data.VEL_MASKED=[];
        data.TEMP=[];
        data.MELTING_LAYER=[];
        %data.SNR=[];

        % Check if LDR_MASKED is available
        try
            velTest=ncread(fileList{1},'LDR_MASKED');
            data.LDR_MASKED=[];
        catch
            data.LDR=[];
        end

        % Load data
        data=read_HCR(fileList,data,startTime,endTime);

        if isfield(data,'LDR_MASKED')
            data.LDR=data.LDR_MASKED;
            data=rmfield(data,'LDR_MASKED');
        end

        ylimits=[0 (max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5];
        plotIn.ylimits=ylimits;

        tempOrig=data.TEMP;

        %% Calculate velocity texture

        pixRadVEL=50;
        velBase=-20;

        data.VEL_MASKED(:,data.elevation<0)=-data.VEL_MASKED(:,data.elevation<0);
        data.VELTEXT=f_velTexture(data.VEL_MASKED,pixRadVEL,velBase);
        
        %% Mask LDR

        data.LDR(isnan(data.DBZ_MASKED))=nan;

        %% Pre process

        disp('Pre processing ...');
        data=preProcessPID(data,convThresh);

        %% Calculate PID

        % HCR
        disp('Getting PID ...');

        plotIn.figdir=[figdir,'debugPlots/'];

        pid_hcr=calc_pid(data,plotIn);

        %% Set areas above melting layer with no LDR to cloud or precip

        data.TEMP=tempOrig;

        smallInds=find((data.MELTING_LAYER>15 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
            (data.DBZ_MASKED<=5 | data.VEL_MASKED<=1));
        pid_hcr(smallInds)=11;

        largeInds=find((data.MELTING_LAYER>15 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
            (data.DBZ_MASKED>5 & data.VEL_MASKED>1));
        pid_hcr(largeInds)=10;

        %% Set low DBZ to cloud liquid

        pid_hcr(pid_hcr>9 & data.DBZ_MASKED<=-30 & data.TEMP>-40)=3;

        %% Set Melting to melting
        pid_hcr(data.MELTING_LAYER==11 | data.MELTING_LAYER==19)=4;

        %% Add supercooled

        disp('Adding supercooled ...')
        pid_hcr=addSupercooled(pid_hcr,data);

        %% Post process

        if postProcess
            pid_hcr=postProcessPID(pid_hcr,data);
        end

        %% Filter

        if whichFilter==1
            pid_hcr=modeFilter(pid_hcr,7,0.7);
        elseif whichFilter==2
            pid_hcr=coherenceFilter(pid_hcr,7,0.7);
        end

        %% Scales and units
        cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];

        units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
            'Melting','Large Frozen','Small Frozen','Precip','Cloud'};
        
        %% Plot PIDs
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        disp('Plotting PID');

        close all
        
        f1=figure('DefaultAxesFontSize',12,'Position',[0 300 2300 1200],'visible','off');
        
        s1=subplot(4,2,1);
        surf(data.time,data.asl./1000,data.DBZ_MASKED,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-40 20]);
        colormap(s1,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR reflectivity (dBZ)']);
        grid on
        
        s3=subplot(4,2,3);
        surf(data.time,data.asl./1000,data.VEL_MASKED,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-5 5]);
        colormap(s3,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR radial velocity (m s^{-1})']);
        grid on
        
        s8=subplot(4,2,5);
        surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-30 -10]);
        colormap(s8,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR linear depolarization ratio (dB)']);
        grid on

        s2=subplot(4,2,2);
        plotMelt=data.MELTING_LAYER;
        meltPlot=nan(size(plotMelt));
        meltPlot(plotMelt==9)=3;
        meltPlot(plotMelt==11)=2;
        meltPlot(plotMelt==19)=1;
        meltPlot(plotMelt==21)=0;
        surf(data.time,data.asl./1000,meltPlot,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        s2.Colormap=[0,1,1;0.5,0.5,0.5;0,0,0;1,0,1];
        caxis([-0.5,3.5]);
        cb=colorbar;
        cb.Ticks=[0,1,2,3];
        cb.TickLabels={'Cold','Melting cold','Melting warm','Warm'};
        ylabel('Altitude (km)');
        title(['Melting Layer']);
        grid on

        s4=subplot(4,2,4);
        jetIn=jet;
        jetTemp=cat(1,jetIn(1:size(jetIn,1)/2,:),repmat([0 0 0],3,1),...
            jetIn(size(jetIn,1)/2+1:end,:));
        surf(data.time,data.asl./1000,data.TEMP,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-15 15]);
        colormap(s4,jetTemp);
        colorbar;
        ylabel('Altitude (km)');
        title(['Temperature (C)']);
        grid on
        
        s6=subplot(4,2,6);
        surf(data.time,data.asl./1000,data.VELTEXT,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([0 5]);
        colormap(s6,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR velocity texture']);
        grid on

        s8=subplot(4,2,8);
        surf(data.time,data.asl./1000,pid_hcr,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 11.5]);
        colormap(s8,cscale_hcr);
        cb=colorbar;
        cb.Ticks=1:11;
        cb.TickLabels=units_str_hcr;
        ylabel('Altitude (km)');
        title(['HCR particle ID']);
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_pid_',...
            datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        
    end
end