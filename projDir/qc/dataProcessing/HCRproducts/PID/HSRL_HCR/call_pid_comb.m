% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.2';
freqData='combined'; % 10hz, 100hz, 2hz, or combined

plotIn.plotMR=0;
plotIn.plotMax=0;

convThresh=4;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

indir=HCRdir(project,quality,qcVersion,freqData);

figdir=[indir(1:end-4),'pidPlotsComb/cases/'];

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
        
        % HCR
        data.HCR_DBZ=[];
        data.HCR_VEL=[];
        data.HCR_LDR=[];
        data.HCR_MELTING_LAYER=[];

        % HSRL
        data.HSRL_Aerosol_Backscatter_Coefficient=[];
        data.HSRL_Particle_Linear_Depolarization_Ratio=[];

        % TEMP
        data.TEMP=[];
                       
        dataVars=fieldnames(data);
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
       
        ylimits=[0 (max(data.asl(~isnan(data.HCR_DBZ)))./1000)+0.5];
        plotIn.ylimits=ylimits;

        tempOrig=data.TEMP;
        %% Calculate velocity texture

        pixRadVEL=10;
        velBase=-20;

        data.HCR_VEL(:,data.elevation<0)=-data.HCR_VEL(:,data.elevation<0);
        data.VELTEXT=f_velTexture(data.HCR_VEL,pixRadVEL,velBase);

        %% Mask LDR and HSRL

        data.HCR_LDR(isnan(data.HCR_DBZ))=nan;
        data.HSRL_Aerosol_Backscatter_Coefficient(isnan(data.HCR_DBZ))=nan;
        data.HSRL_Particle_Linear_Depolarization_Ratio(isnan(data.HCR_DBZ))=nan;

        %% Pre process

        disp('Pre processing ...');
        data=preProcessPIDcomb(data,convThresh);

        %% Calculate PID

        % HCR
        disp('Getting PID ...');

        plotIn.figdir=[figdir,'debugPlots/'];

        pid_hcr_hsrl=calc_pid_comb(data,plotIn);

        %% Set areas above melting layer with no LDR to cloud or precip

        data.TEMP=tempOrig;

        smallInds=find((data.HCR_MELTING_LAYER>15 | (isnan(data.HCR_MELTING_LAYER) & data.TEMP<0)) & isnan(data.HCR_LDR) & isnan(data.HSRL_Particle_Linear_Depolarization_Ratio) & ...
            (data.HCR_DBZ<=5 | data.HCR_VEL<=1));
        pid_hcr_hsrl(smallInds)=11;

        largeInds=find((data.HCR_MELTING_LAYER>15 | (isnan(data.HCR_MELTING_LAYER) & data.TEMP<0)) & isnan(data.HCR_LDR) & isnan(data.HSRL_Particle_Linear_Depolarization_Ratio) & ...
            (data.HCR_DBZ>5));
        pid_hcr_hsrl(largeInds)=10;

        %% Set low DBZ to cloud liquid

        pid_hcr_hsrl(pid_hcr_hsrl>9 & data.HCR_DBZ<=-30)=3;

        %% Set Melting to melting
        pid_hcr_hsrl(data.HCR_MELTING_LAYER==11 | data.HCR_MELTING_LAYER==19)=4;

        %% Add supercooled

        disp('Adding supercooled ...')
        pid_hcr_hsrl=addSupercooledComb(pid_hcr_hsrl,data);

        %% Post process

        if postProcess
            disp('Post processing ...');
            pid_hcr_hsrl=postProcessPIDcomb(pid_hcr_hsrl,data);
        end

        %% Filter

        if whichFilter==1
            pid_hcr_hsrl=modeFilter(pid_hcr_hsrl,7,0.7);
        elseif whichFilter==2
            pid_hcr_hsrl=coherenceFilter(pid_hcr_hsrl,7,0.7);
        end

        %% Scales and units
        cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1; 0,0,0; 0.5,0.5,0.5];
        
        units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid',...
            'Mixed Phase','Large Frozen','Small Frozen','Precip','Cloud'};
        
        %% Plot PIDs
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        disp('Plotting PID');

        close all
        
        f1=figure('DefaultAxesFontSize',12,'Position',[0 300 2300 1200],'visible','off');
        
        s1=subplot(5,2,1);
        surf(data.time,data.asl./1000,data.HCR_DBZ,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-40 20]);
        colormap(s1,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR reflectivity (dBZ)']);
        grid on
        
        s3=subplot(5,2,3);
        surf(data.time,data.asl./1000,data.HCR_VEL,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-5 5]);
        colormap(s3,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR radial velocity (m s^{-1})']);
        grid on
        
        s5=subplot(5,2,5);
        surf(data.time,data.asl./1000,data.HCR_LDR,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-30 -10]);
        colormap(s5,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR linear depolarization ratio (dB)']);
        grid on

        s7=subplot(5,2,7);
        surf(data.time,data.asl./1000,log10(data.HSRL_Aerosol_Backscatter_Coefficient),'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        %caxis([-30 -10]);
        colormap(s7,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HSRL aerosol backscatter coefficient (m^{-1} sr^{-1})']);
        grid on
                
        s2=subplot(5,2,2);
        plotMelt=data.HCR_MELTING_LAYER;
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
        
        s4=subplot(5,2,4);
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
        
        s6=subplot(5,2,6);
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

        s8=subplot(5,2,8);
        surf(data.time,data.asl./1000,data.HSRL_Particle_Linear_Depolarization_Ratio,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([0 0.6]);
        colormap(s8,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HSRL particle linear depolarization ratio']);
        grid on

        s10=subplot(5,2,10);
        surf(data.time,data.asl./1000,pid_hcr_hsrl,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 11.5]);
        colormap(s10,cscale_hcr);
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