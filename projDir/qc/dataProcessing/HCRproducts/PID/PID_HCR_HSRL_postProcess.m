% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='combined'; % 10hz, 100hz, 2hz, or combined

ylimits=[0 3];

plotComp=1; % 1 to plot comparison plot of HCR vs HSRL
plotFields=1; % 1 to plot input fields
whichFilter=1; % 0: no filter, 1: mode filter, 2: coherence filter

figdir='/home/romatsch/plots/HCR/pid/postProcess/';

%indir=HCRdir(project,quality,freqData);
%indir=HCRdirWFH(project,quality,freqData);
indir='/run/media/romatsch/RSF0006/rsf/meltingLayer/socrates/combined/';

% Loop through cases
casefile=['~/git/HCR_configuration/projDir/qc/dataProcessing/HCRproducts/caseFiles/pid_',project,'.txt'];

caseList=readtable(casefile);
caseStart=datetime(caseList.Var1,caseList.Var2,caseList.Var3, ...
    caseList.Var4,caseList.Var5,0);
caseEnd=datetime(caseList.Var6,caseList.Var7,caseList.Var8, ...
    caseList.Var9,caseList.Var10,0);

for aa=3:length(caseStart)
    
    disp(['Case ',num2str(aa),' of ',num2str(length(caseStart))]);
    
    startTime=caseStart(aa);
    endTime=caseEnd(aa);
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        %% Load data
        
        data=[];
        
        %HCR data
        data.HCR_DBZ=[];
        data.HCR_VEL=[];
        data.HCR_WIDTH=[];
        data.HCR_LDR=[];
        data.TEMP=[];
        data.MELTING_LAYER=[];
        
        %HSRL data
        data.HSRL_Aerosol_Backscatter_Coefficient=[];
        data.HSRL_Volume_Depolarization=[];
        %data.HSRL_Aerosol_Extinction_Coefficient=[];
        
        dataVars=fieldnames(data);
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
        
        data.asl=data.asl./1000;
        data.temp=data.TEMP+273.15;
        
        %% Calculate HSRL PID
        
        % HSRL
        backscatLog = real(log10(data.HSRL_Aerosol_Backscatter_Coefficient));
        %extLog = real(log10(data.HSRL_Aerosol_Extinction_Coefficient));
        depolLog = real(log10(data.HSRL_Volume_Depolarization));
        %lidarRatio=10.^(extLog-backscatLog);
        vol_depol=data.HSRL_Volume_Depolarization./(2-data.HSRL_Volume_Depolarization);
        lin_depol=vol_depol./(2-vol_depol);
        
        pid_hsrl=calc_pid_hsrl_postProcess(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
        pid_hsrl(isnan(data.HSRL_Aerosol_Backscatter_Coefficient))=nan;
        
        if whichFilter==1
            pid_hsrl=modeFilter(pid_hsrl,7,0.7);
        elseif whichFilter==2
            pid_hsrl=coherenceFilter(pid_hsrl,7,0.7);
        end        
          
        %% Calculate HCR without attenuation correction
        
        [pid_hcr]=calc_pid_hcr_postProcess(data.HCR_DBZ,data);
        pid_hcr(isnan(data.HCR_DBZ))=nan;
             
        % Combined from merging hcr and hsrl pid
        pid_comb=combine_pid_hcr_hsrl_postProcess(pid_hcr,pid_hsrl);
        
%         % Combined by using both data sets in one process
%         pid_comb2=calc_pid_direct_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
%             data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
        
        %% Calculate attenuation correction
        
        Z_95_lin=10.^(data.HCR_DBZ*0.1);
        Z_95_lin(data.HCR_DBZ < -200)=0.;
        
        % Mask out non liquid data
        liqMeltInds=find(pid_comb==1 | pid_comb==2 | pid_comb==3 | pid_comb==4);
        Z_95_lin(liqMeltInds)=nan;
        
        wt_coef=nan(size(data.HCR_DBZ));
        wt_exp=nan(size(data.HCR_DBZ));
        
        wt_coef(data.HCR_DBZ < - 20)=20.;
        wt_exp(data.HCR_DBZ < - 20)=0.52;
        wt_coef(-20 <data.HCR_DBZ <-15 )=1.73;
        wt_exp(-20 <data.HCR_DBZ < -15 )=0.15;
        wt_coef(data.HCR_DBZ > -15)=0.22;
        wt_exp(data.HCR_DBZ > -15)=0.68;
        
        att_cumul=2.*0.0192*cumsum((wt_coef.*Z_95_lin.^wt_exp),1,'omitnan');
        att_cumul(data.HCR_DBZ < -200)=NaN;
        dBZ_cor_all=data.HCR_DBZ+att_cumul;
        
        % Replace dBZ values with attenuation corrected values in liquid and
        % melting regions
        dBZ_cor=data.HCR_DBZ;
        dBZ_cor(liqMeltInds)=dBZ_cor_all(liqMeltInds);
        dBZ_cor(isnan(data.HCR_DBZ))=nan;
        
        %% Calculate PID with attenuation correction
        
        % HCR
        [pid_hcr_cor]=calc_pid_hcr_postProcess(dBZ_cor,data);
        pid_hcr_cor(isnan(dBZ_cor))=nan;
        
        if whichFilter==1
            pid_hcr_cor=modeFilter(pid_hcr_cor,7,0.7);
        elseif whichFilter==2
            pid_hcr_cor=coherenceFilter(pid_hcr_cor,7,0.7);
        end
        
        % Combined from merging hcr and hsrl pid
        [pid_comb_cor which_pid]=combine_pid_hcr_hsrl_postProcess(pid_hcr_cor,pid_hsrl);
        
        which_pid(isnan(pid_comb_cor))=nan;
        disagree_pid=zeros(size(pid_comb_cor));
        disagree_pid(find(pid_hcr_cor~=pid_hsrl))=1;
        disagree_pid(find(~isnan(pid_hcr_cor) & isnan(pid_hsrl)))=2;
        disagree_pid(find(isnan(pid_hcr_cor) & ~isnan(pid_hsrl)))=3;
        disagree_pid(isnan(pid_comb_cor))=nan;
        
        % Combined by using both data sets in one process
        %         pid_comb2_cor=calc_pid_direct_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
        %             dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
        
        %% Scales and units
        cscale_hsrl=[0,0,1.0;0,1,0;1,0.67,0;1,0,1;0,1,1;1,0.67,0];
        cscale_hcr=[0,0,1.0; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
        cscale_comb=[0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
        
        units_str_hsrl={'Cloud liquid','Drizzle',...
            'Aerosol1','SLW','Ice crystals','Aerosol2'};
        units_str_hcr={'Cloud liquid','Drizzle',...
            'Rain','SLW','Ice crystals','Snow','Wet snow/rimed ice'};
        units_str_comb={'Cloud liquid','Drizzle','Rain',...
            'SLW','Ice crystals','Snow','Wet snow/rimed ice','Aerosols'};
        
        %% Plot PIDs
        
        timeMat=repmat(data.time,size(data.TEMP,1),1);
        
        close all
        f1=figure('DefaultAxesFontSize',12,'Position',[400 300 1200 900]);
        
        s1=subplot(3,1,1);
        surf(data.time,data.asl,pid_hsrl,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 6.5]);
        colormap(s1,cscale_hsrl);
        cb=colorbar;
        cb.Ticks=1:6;
        cb.TickLabels=units_str_hsrl;
        ylabel('Altitude (km)');
        title(['HSRL particle ID']);
        
        s2=subplot(3,1,2);
        surf(data.time,data.asl,pid_hcr_cor,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 7.5]);
        colormap(s2,cscale_hcr);
        cb=colorbar;
        cb.Ticks=1:7;
        cb.TickLabels=units_str_hcr;
        ylabel('Altitude (km)');
        title(['HCR particle ID']);
        
        s3=subplot(3,1,3);
        surf(data.time,data.asl,pid_comb_cor,'edgecolor','none');
        view(2);
        ylim(ylimits);
        xlim([data.time(1),data.time(end)]);
        caxis([.5 8.5]);
        colormap(s3,cscale_comb);
        cb=colorbar;
        cb.Ticks=1:8;
        cb.TickLabels=units_str_comb;
        ylabel('Altitude (km)');
        title(['Combined particle ID']);
        
        set(gcf,'PaperPositionMode','auto')
        print(f1,[figdir,project,'_pid_',...
            datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_pid'],'-dpng','-r0')
        
        %% Plot comparison plot
        if plotComp
            f2=figure('DefaultAxesFontSize',12,'Position',[400 300 1200 900]);
            
            s1=subplot(3,1,1);
            fig1=surf(data.time,data.asl,pid_comb_cor,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([.5 8.5]);
            colormap(s1,cscale_comb);
            cb=colorbar;
            cb.Ticks=1:8;
            cb.TickLabels=units_str_comb;
            ylabel('Altitude (km)');
            title(['Combined particle ID']);
                        
            s2=subplot(3,1,2);
            fig2=surf(data.time,data.asl,which_pid,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([0 1]);
            colormap(s2,[0 0 1; 1 0 0]);
            cb2=colorbar;
            cb2.Ticks=[0.25,0.75];
            cb2.TickLabels={'HCR','HSRL'};
            ylabel('Altitude (km)');
            title(['PID source']);
            
            s3=subplot(3,1,3);
            fig3=surf(data.time,data.asl,disagree_pid,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([0 4]);
            colormap(s3,[0 1 1; 1 1 0;0 0 1; 1 0 0]);
            cb3=colorbar;
            cb3.Ticks=[0.5,1.5,2.5,3.5];
            cb3.TickLabels={'Same','Different','HCR only','HSRL only'};
            ylabel('Altitude (km)');
            title(['PID agreement']);
            
            set(gcf,'PaperPositionMode','auto')
            print(f2,[figdir,project,'_pid_',...
                datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_debug'],'-dpng','-r0')
            
        end
        
        %% Plot input fields
        if plotFields
            f3=figure('DefaultAxesFontSize',12,'Position',[0 300 1700 900]);
            
            s1=subplot(4,2,1);
            surf(data.time,data.asl,data.HCR_DBZ,'edgecolor','none');
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
            surf(data.time,data.asl,data.HCR_VEL,'edgecolor','none');
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
            surf(data.time,data.asl,data.HCR_LDR,'edgecolor','none');
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
            surf(data.time,data.asl,data.HCR_WIDTH,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([0 2]);
            colormap(s7,jet);
            colorbar;
            ylabel('Altitude (km)');
            title(['HCR spectrum width (m s^{-1})']);
            grid on
            
            s2=subplot(4,2,2);
            plotMelt=data.MELTING_LAYER;
            plotMelt(~isnan(plotMelt) & plotMelt<20)=10;
            plotMelt(~isnan(plotMelt) & plotMelt>=20)=20;
            surf(data.time,data.asl,plotMelt,'edgecolor','none');
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
            surf(data.time,data.asl,data.TEMP,'edgecolor','none');
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
            fig3=surf(data.time,data.asl,backscatLog,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([-7 -2]);
            colormap(s6,jet);
            colorbar;
            ylabel('Altitude (km)');
            title(['HSRL log10 aerosol backscatter coefficient (m^{-1} sr^{-1})']);
            
            s8=subplot(4,2,8);
            fig3=surf(data.time,data.asl,data.HSRL_Volume_Depolarization,'edgecolor','none');
            view(2);
            ylim(ylimits);
            xlim([data.time(1),data.time(end)]);
            caxis([0 0.8]);
            colormap(s8,jet);
            colorbar;
            ylabel('Altitude (km)');
            title(['HSRL volume depolarization']);
            
            set(gcf,'PaperPositionMode','auto')
            print(f3,[figdir,project,'_pid_',...
                datestr(data.time(1),'yyyymmdd_HHMMSS'),'_to_',datestr(data.time(end),'yyyymmdd_HHMMSS'),'_fields'],'-dpng','-r0')
        end
        
    end
end