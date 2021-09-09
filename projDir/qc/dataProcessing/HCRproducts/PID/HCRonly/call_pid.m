% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
qcVersion='v2.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined

ylimits=[0 15];
plotIn.plotMR=0;
plotIn.plotMax=0;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

if strcmp(project,'otrec')
    indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
elseif strcmp(project,'socrates')
    indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
end

figdir=[indir(1:end-5),'pidPlots/cases/'];

plotIn.figdir=figdir;
plotIn.ylimits=ylimits;

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
        data.DBZ=[];
        data.VEL_CORR=[];
        data.WIDTH=[];
        data.LDR=[];
        data.TEMP=[];
        data.MELTING_LAYER=[];
        data.FLAG=[];
        
        dataVars=fieldnames(data);
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
        
        % Mask with FLAG
        data.DBZ(data.FLAG>1)=nan;
        data.VEL_CORR(data.FLAG>1)=nan;
        data.WIDTH(data.FLAG>1)=nan;
        data.LDR(data.FLAG>1)=nan;
        data.TEMP(data.FLAG>1)=nan;
        data.MELTING_LAYER(data.FLAG>1)=nan;
        
        %% Correct for attenuation
        
        Z_lin=10.^(data.DBZ*0.1);
        
        % Mask out non liquid data
        liqMeltInds=find(data.MELTING_LAYER<20);
        Z_lin(~liqMeltInds)=nan;
        
        wt_coef=nan(size(data.DBZ));
        wt_exp=nan(size(data.DBZ));
        
        wt_coef(data.DBZ < - 20)=20.;
        wt_exp(data.DBZ < - 20)=0.52;
        wt_coef(-20 <data.DBZ <-15 )=1.73;
        wt_exp(-20 <data.DBZ < -15 )=0.15;
        wt_coef(data.DBZ > -15)=0.22;
        wt_exp(data.DBZ > -15)=0.68;
        
        att_cumul=2.*0.0192*cumsum((wt_coef.*Z_lin.^wt_exp),1,'omitnan');
        dBZ_cor_all=data.DBZ+att_cumul;
        
        % Replace dBZ values with attenuation corrected values in liquid and
        % melting regions
        dBZ_cor=data.DBZ;
        dBZ_cor(liqMeltInds)=dBZ_cor_all(liqMeltInds);
        
        %% Calculate PID with attenuation correction
        
        % HCR
        disp('Getting PID');
        [pid_hcr]=calc_pid(dBZ_cor,data,postProcess,plotIn);
        
        if whichFilter==1
            pid_hcr=modeFilter(pid_hcr,7,0.7);
        elseif whichFilter==2
            pid_hcr=coherenceFilter(pid_hcr,7,0.7);
        end
        
        
        %% Scales and units
        cscale_hcr=[1,0,0; 1,0.6,0.47; 0,1,0; 0,0.7,0; 0,0,1; 1,0,1; 0.5,0,0; 1,1,0; 0,1,1];
        
        units_str_hcr={'Rain','Supercooled Rain','Drizzle','Supercooled Drizzle','Cloud Liquid','Supercooled Cloud Liquid','Mixed Phase','Large Frozen','Small Frozen'};
        
        %% Plot PIDs
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        disp('Plotting PID');
        
        close all
        
        f1=figure('DefaultAxesFontSize',12,'Position',[0 300 2300 1200],'visible','off');
        
        s1=subplot(4,2,1);
        surf(data.time,data.asl./1000,data.DBZ,'edgecolor','none');
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
        surf(data.time,data.asl./1000,data.VEL_CORR,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-5 5]);
        colormap(s3,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR radial velocity (m s^{-1})']);
        grid on
        
        s5=subplot(4,2,5);
        surf(data.time,data.asl./1000,data.LDR,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([-30 -5]);
        colormap(s5,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR linear depolarization ratio (dB)']);
        grid on
        
        s7=subplot(4,2,7);
        surf(data.time,data.asl./1000,data.WIDTH,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([0 1]);
        colormap(s7,jet);
        colorbar;
        ylabel('Altitude (km)');
        title(['HCR spectrum width (m s^{-1})']);
        grid on
        
        s2=subplot(4,2,2);
        plotMelt=data.MELTING_LAYER;
        plotMelt(~isnan(plotMelt) & plotMelt<20)=10;
        plotMelt(~isnan(plotMelt) & plotMelt>=20)=20;
        surf(data.time,data.asl./1000,plotMelt,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        %caxis([0 2]);
        colormap(s2,[1 0 1;1 1 0]);
        c=colorbar;
        c.Visible='off';
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
        
        s5=subplot(4,2,6);
        surf(data.time,data.asl./1000,pid_hcr,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 9.5]);
        colormap(s5,cscale_hcr);
        cb=colorbar;
        cb.Ticks=1:9;
        cb.TickLabels=units_str_hcr;
        ylabel('Altitude (km)');
        title(['HCR particle ID']);
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_pid_',...
            datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS')],'-dpng','-r0')
        
        
    end
end