% Read HCR HSRL combined data

clear all
close all

% startTime=datetime(2018,1,24,3,59,00);
% endTime=datetime(2018,1,24,4,00,00);

% startTime=datetime(2018,1,23,00,00,0);
% endTime=datetime(2018,1,23,01,00,0);

% startTime=datetime(2018,1,24,3,50,0); %BAMS Jeff Stith
% %  startTime=datetime(2018,1,24,4,01,0); %BAMS Jeff Stith
%  endTime=datetime(2018,1,24,4,05,0); %BAMS Jeff Stith
 
 
 startTime=datetime(2018,1,24,01,09,30); %Wang_Rauber
 endTime=datetime(2018,1,24,01,12,30); %Wang_Rauber
%
%
% startTime=datetime(2015,7,24,19,15,0);
% endTime=datetime(2015,7,24,19,20,0);

%  startTime=datetime(2018,2,20,3,19,0);% JGR
%  endTime=datetime(2018,2,20,3,24,0); %  JGR

ylimits=[3.0 6.0];%ylimits=[-0.0 2.5];

plotlidars=1; % 1 to plot lidar data, 0 to not plot lidar
plotradars=1; % 1 to plot radar data, 0 to not plot radar
filltemp=1; % 1 to fill in temperature data with workaround to masking issue

indir='/Volumes/RSF-Vivek/SOCRATES/HCR_HSRL_qc2_RF04_20180123_230524_to_20180124_060037/';
%indir='/Volumes/RSF-Vivek/CSET/2hz_merge_15min.cset/';
formatOut = 'yyyymmdd_HHMM';

fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

if ~isempty(fileList)
    
    %HCR data
    data.HCR_DBZ=[];
    data.HCR_VEL=[];
    data.HCR_WIDTH=[];
    data.HCR_LDR=[];
    data.TEMP=[];
    
    %HSRL data
    data.HSRL_Aerosol_Backscatter_Coefficient=[];
    data.HSRL_Volume_Depolarization=[];
    data.HSRL_Aerosol_Extinction_Coefficient=[];   
    
    dataVars=fieldnames(data);   
    
    % Load data
    data=read_HCR(fileList,data,startTime,endTime);
    
    % Check if all variables were found
    for ii=1:length(dataVars)
        if ~isfield(data,dataVars{ii})
            dataVars{ii}=[];
        end
    end
       
    Z_95_lin=10.^(data.HCR_DBZ*0.1);
    
    
    wt_coef=zeros(size(data.HCR_DBZ));
    wt_exp=zeros(size(data.HCR_DBZ));
    Z_95_lin(data.HCR_DBZ < -200)=0.;
    DBZ_temp=data.HCR_DBZ;
    wt_coef(data.HCR_DBZ < - 20)=20.;
    wt_exp(data.HCR_DBZ < - 20)=0.52;
    wt_coef(-20 <data.HCR_DBZ <-15 )=1.73;
    wt_exp(-20 <data.HCR_DBZ < -15 )=0.15;
    wt_coef(data.HCR_DBZ > -15)=0.22;
    wt_exp(data.HCR_DBZ > -15)=0.68;
    att_cumul=2.*0.0192*cumsum((wt_coef.*Z_95_lin.^wt_exp),2,'omitnan');
    att_cumul(data.HCR_DBZ < -200)=NaN;
    dBZ_cor=data.HCR_DBZ+att_cumul;
    Z_95_lin_cor=10.^(dBZ_cor*0.1);
    
    %% Calculate variables
    
    data.temp=data.TEMP+273.15;
    
    close all
    %calculate above sea level altitudes
    asl=nan(size(data.range));
    downInd=find(data.elevation<0);
    upInd=find(data.elevation>=0);
    asl(:,downInd)=-1*((data.range(:,downInd).*cosd(abs(data.elevation(downInd))-90)./1000)-data.altitude(downInd)./1000);
    asl(:,upInd)=data.range(:,upInd).*cosd(abs(data.elevation(upInd))-90)./1000+data.altitude(upInd)./1000;
    
    backscatLog = real(log10(data.HSRL_Aerosol_Backscatter_Coefficient));
    extLog = real(log10(data.HSRL_Aerosol_Extinction_Coefficient));
    depolLog = real(log10(data.HSRL_Volume_Depolarization));
    
    lidarRatio=10.^(extLog-backscatLog);
    
    %Particle ID
    vol_depol=data.HSRL_Volume_Depolarization./(2-data.HSRL_Volume_Depolarization);
    lin_depol=vol_depol./(2-vol_depol);
    %     pid=aer_cld2c(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
    pid=aer_cld2c_ver5(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
    
    %v[pidhcr,m]=precip_cld_a(data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL_CORR,data.HCR_WIDTH_CORR,data.temp);
    
    %[pidhcr,m]=precip_cld_a(data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    
    [pidhcr,m]=precip_cld_a(dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    
    %vpidhcr(isnan(data.HCR_DBZ))=nan;
    pid(isnan(data.HSRL_Aerosol_Backscatter_Coefficient))=nan;
    pidhcr(isnan(dBZ_cor))=nan;
    
    pid_comb=combine_pid_hcr_hsrl(pidhcr,pid,data.temp);
    
%     pid_comb2=precip_aerosol(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
%         data.HCR_DBZ,data.HCR_LDR, data.HCR_VEL,data.HCR_WIDTH,data.temp);
    
    pid_comb2=precip_aerosol(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
        dBZ_cor,data.HCR_LDR, data.HCR_VEL,data.HCR_WIDTH,data.temp);
    %% Plot
    close all
    if plotlidars==1
        f1=figure('DefaultAxesFontSize',12);
        set(f1,'Position',[400 300 1000 700]);
        
        colormap(jet);
        
        subplot(4,1,1)
        fig1=surf(data.time,asl,backscatLog);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([-7 -2]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title(['Log10 aerosol backscatter coefficient (m^-1 sr^-1)'],'interpreter','none');
        
        
        
        %         subplot(4,1,2)
        %         hold on
        %         fig1=surf(data.time,asl,extLog);
        %         fig1.EdgeColor='none';
        %         ylim(ylimits);
        %         xlim([startTime,endTime]);
        %         caxis([-6 -1]);
        %         colorbar;
        %         view(2);
        %         ylabel('Altitude (km)');
        %         title(['Log10 aerosol extinction coefficient (m^-1)'],'interpreter','none');
        
        
        subplot(4,1,2)
        fig1=surf(data.time,asl,data.HSRL_Volume_Depolarization);
        fig1.EdgeColor='none';
        ylim(ylimits);
        caxis([0 0.8]);
        xlim([startTime,endTime]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title(['Aerosol Linear depolarization ratio'],'interpreter','none');
        
        subplot(4,1,3)
        fig1=surf(data.time,asl,data.temp);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([260 280]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title(['Temp, K'],'interpreter','none');
   %                1 .     2 .      3 .       4 .    5 .    6 .     7
        %cscale=[1,1,1; 1,0.67,0; 1,0.67,0; 0,1,1; 1,0,1; 0,0,0.5; 0,1,0];
        cscale=[1,1,1;0,0,0.5;0,1,0;1,0.67,0;1,0,1;0,1,1;1,0.67,0];
%         units_str={'No Signal', 'Aerosol','', 'Ice', 'SLD', 'CLD', 'Drizzle'};
       % units_str={'No Signal', 'Aerosol','', 'Ice', 'Droplets below freezing', 'Droplets above freezing', 'Drizzle'};
        units_str={'No Signal', 'Droplets above freezing','Drizzle', 'Aerosol1', 'Droplets below freezing', 'Ice crystals', 'Aerosol2'};
        subplot(4,1,4)
        fig1=surf(data.time,asl,pid);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([.5 7.5]);
        cb=colorbar;
        colormap(gca,cscale);
        cb.Ticks=1:7;
        cb.TickLabels=units_str;
        view(2);
        ylabel('Altitude (km)');
        title(['Particle ID HSRL'],'interpreter','none');
        
    end
    
    if plotradars==1
        
        f2=figure('DefaultAxesFontSize',16);
        set(f2,'Position',[600 270 1000 700]);
        set(gca,'fontname','Timesnewroman','FontSize',14','fontweight','b');
        %         colormap(jet);
        
        % Vel Raw
        subplot(4,1,1)
        set(gca,'FontSize',18);
        set(gca,'FontWeight','bold');
        fig1=surf(data.time,asl,data.HCR_DBZ);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]); set(gca,'xticklabel',[]); %set(gca,'xtick',[])
        caxis([-30 15]);
        colorbar;
        view(2);
        ylabel('Altitude (km)','fontsize',10,'fontweight','b');
        title('Radar reflectivity (dBZ)  ', ...
            'interpreter','none');
        colormap(gca,cid_cmap2)
        set(gca,'TickLabelInterpreter','none');
        set(gca,'fontweight','bold','fontsize',16);
        
        %         subplot(5,1,2)
        %         %hold on
        %         fig1=surf(data.time,asl,data.HCR_LDR);
        %         fig1.EdgeColor='none';
        %         ylim(ylimits);
        %         xlim([startTime,endTime]);set(gca,'xticklabel',[]);
        %         caxis([-25 -15]);
        %         colorbar;
        %         view(2);
        %         ylabel('Altitude (km)','fontsize',10,'fontweight','b');
        %         title(['Radar linear depolarization ratio (dB)'],'interpreter','none');
        %         colormap(gca,cid_cmap2);
        %         set(gca,'TickLabelInterpreter','none');
        %         set(gca,'fontweight','bold','fontsize',12);
        
        
        
        %         subplot(4,1,3)
        %         hold on
        %         fig1=surf(data.time,asl,data.HCR_VEL_CORR);
        %         fig1.EdgeColor='none';
        %         ylim(ylimits);
        %         xlim([startTime,endTime]);
        %         caxis([-2 6]);
        %         colorbar;
        %         view(2);
        %         ylabel('Altitude (km)');
        %         title(['Vel m/s)'],'interpreter','none');
        %         colormap(gca,cid_cmap)
        
        subplot(4,1,2)
        % hold on
        fig1=surf(data.time,asl,backscatLog);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);set(gca,'xticklabel',[]);
        caxis([-7 -2]);
        colorbar;
        view(2);
        ylabel('Altitude (km)','fontsize',10,'fontweight','b');
        title(['Lidar backscatter log_{10}[{\beta}] m^{-1} Sr^{-1}']);
        set(gca,'TickLabelInterpreter','none');
        set(gca,'fontweight','bold','fontsize',16);
        %         subplot(4,1,4)
        %         hold on
        %         fig1=surf(data.time,asl,data.HCR_WIDTH_CORR);
        %         fig1.EdgeColor='none';
        %         ylim(ylimits);
        %         xlim([startTime,endTime]);
        %         caxis([0 2]);
        %         colorbar;
        %         view(2);
        %         ylabel('Altitude (km)');
        %         title(['Spectrum width m/s)'],'interpreter','none');
        %         colormap(gca,cid_cmap2)
        
        subplot(4,1,3)
        %hold on
        fig1=surf(data.time,asl,data.HSRL_Volume_Depolarization);
        fig1.EdgeColor='none';
        ylim(ylimits);
        caxis([0 0.8]);
        xlim([startTime,endTime]);set(gca,'xticklabel',[]);
        colorbar;
        view(2);
        ylabel('Altitude (km)','fontsize',10,'fontweight','b');
        title(['Lidar depolarization ratio'],'interpreter','none');
        set(gca,'TickLabelInterpreter','none');
        set(gca,'fontweight','bold','fontsize',16);
        
        subplot(4,1,4)
        
        
        cscale=[1,1,1; 0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
        units_str={'No signal','Droplets above freezing', 'Drizzle','Rain', 'Droplets below freezing',...
            'Ice crystals', 'Rimed ice', 'Wet snow', 'Aerosol'};
        fig1=surf(data.time,asl,pid_comb);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([.5 9.5]);
        cb=colorbar;
        colormap(gca,cscale);
        cb.Ticks=1:9;
        cb.TickLabels=units_str;
        view(2);
        ylabel('Altitude (km)','fontsize',10,'fontweight','b');
        xlabel('Time (UTC)','fontsize',10,'fontweight','b')
        title(['Particle ID Combined'],'interpreter','none');
        set(gca,'TickLabelInterpreter','none');
        set(gca,'fontweight','bold','fontsize',16);
        
        cscale=[1,1,1; 0,0,.5; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
        
        units_str={'No signal','Droplets above freezing', 'Drizzle','Rain', 'Droplets below freezing', 'Ice crystals', 'Snow', 'Wet snow/rimed ice'};
        
%         subplot(5,1,5)
%         
%         
%         cscale=[1,1,1; 0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
%         units_str={'No signal','Droplets above freezing', 'Drizzle','Rain', 'Droplets below freezing',...
%             'Ice crystals', 'Rimed ice', 'Wet snow', 'Aerosol'};
%         fig1=surf(data.time,asl,pid_comb2);
%         fig1.EdgeColor='none';
%         ylim(ylimits);
%         xlim([startTime,endTime]);
%         caxis([.5 9.5]);
%         cb=colorbar;
%         colormap(gca,cscale);
%         cb.Ticks=1:9;
%         cb.TickLabels=units_str;
%         view(2);
%         ylabel('Altitude (km)','fontsize',10,'fontweight','b');
%         xlabel('Time (UTC)','fontsize',10,'fontweight','b')
%         title(['Particle ID Combined'],'interpreter','none');
%         set(gca,'TickLabelInterpreter','none');
%         set(gca,'fontweight','bold','fontsize',16);
        
        cscale=[1,1,1; 0,0,.5; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
        
        units_str={'No signal','Droplets above freezing', 'Drizzle','Rain', 'Droplets below freezing', 'Ice crystals', 'Snow', 'Wet snow/rimed ice'};
        
        f3=figure('DefaultAxesFontSize',16);
        set(f3,'Position',[700 260 1000 700]);
        
        %     colormap(jet);
        
        subplot(5,1,1)
        fig1=surf(data.time,asl,data.HCR_DBZ);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([-30 15]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title('Reflectivity (dBZ)  ', ...
            'interpreter','none');
        colormap(gca,cid_cmap2)
        
        subplot(5,1,2)
        fig1=surf(data.time,asl,data.HCR_VEL);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([-2 4]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title(['Vel m/s)'],'interpreter','none');
        colormap(gca,cid_cmap)
        
        subplot(5,1,3)
        fig1=surf(data.time,asl,data.HCR_LDR);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([-30 -5]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title([' LDR dB)'],'interpreter','none');
        colormap(gca,cid_cmap)
        
        subplot(5,1,4)
        fig1=surf(data.time,asl,data.HCR_WIDTH);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([0 2]);
        colorbar;
        view(2);
        ylabel('Altitude (km)');
        title(['Spectrum width m/s)'],'interpreter','none');
        colormap(gca,cid_cmap)
        
        subplot(5,1,5)
        fig1=surf(data.time,asl,pidhcr);
        %     fig1=surf(data.time,asl,pid_comb);
        fig1.EdgeColor='none';
        ylim(ylimits);
        xlim([startTime,endTime]);
        caxis([.5 8.5]);
        cb=colorbar;
        colormap(gca,cscale);
        cb.Ticks=1:8;
        cb.TickLabels=units_str;
        view(2);
        ylabel('Altitude (km)');
        title(['Particle ID HCR'],'interpreter','none');
        
        
    end
    
    plot_pids;
end