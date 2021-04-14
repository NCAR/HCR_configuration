% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='2hzMerged'; % 10hz, 100hz, or 2hz

startTime=datetime(2018,1,29,1,30,0); %Wang_Rauber
endTime=datetime(2018,1,29,02,0,0); %Wang_Rauber

figdir='/home/romatsch/plots/HCR/pid/noDrizzle/';

ylimits=[0 3];

plotlidars=0; % 1 to plot lidar data, 0 to not plot lidar
plotradars=0; % 1 to plot radar data, 0 to not plot radar

%indir='/Volumes/RSF-Vivek/SOCRATES/HCR_HSRL_qc2_RF04_20180123_230524_to_20180124_060037/';
%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/',project,'/'];

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
    
    data.asl=data.asl./1000;
    data.temp=data.TEMP+273.15;
    
    %% Calculate HSRL PID
    
    % HSRL
    backscatLog = real(log10(data.HSRL_Aerosol_Backscatter_Coefficient));
    extLog = real(log10(data.HSRL_Aerosol_Extinction_Coefficient));
    depolLog = real(log10(data.HSRL_Volume_Depolarization));
    lidarRatio=10.^(extLog-backscatLog);
    vol_depol=data.HSRL_Volume_Depolarization./(2-data.HSRL_Volume_Depolarization);
    lin_depol=vol_depol./(2-vol_depol);
    
    pid_hsrl=calc_pid_hsrl_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,data.temp);
    pid_hsrl(isnan(data.HSRL_Aerosol_Backscatter_Coefficient))=nan;
    pid_hsrl(isnan(pid_hsrl))=1;
    
    %% Calculate HCR without attenuation correction
    
    [pid_hcr]=calc_pid_hcr_clean_eff(data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    pid_hcr(isnan(data.HCR_DBZ))=nan;
    pid_hcr(isnan(pid_hcr))=1;
    
    % Combined from merging hcr and hsrl pid
    pid_comb=combine_pid_hcr_hsrl_clean(pid_hcr,pid_hsrl);
    
    % Combined by using both data sets in one process
    %         pid_comb2=calc_pid_direct_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
    %             data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    %
    %% Calculate attenuation correction
    
    Z_95_lin=10.^(data.HCR_DBZ*0.1);
    Z_95_lin(data.HCR_DBZ < -200)=0.;
    
    % Mask out non liquid data
    liqMeltInds=find(pid_comb==2 | pid_comb==3 | pid_comb==4 | pid_comb==5);
    Z_95_lin(liqMeltInds)=nan;
    
    %DBZ_temp=data.HCR_DBZ;
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
    [pid_hcr_cor]=calc_pid_hcr_clean_eff(dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    pid_hcr_cor(isnan(dBZ_cor))=nan;
    pid_hcr_cor(isnan(pid_hcr_cor))=1;
    
    % Combined from merging hcr and hsrl pid
    pid_comb_cor=combine_pid_hcr_hsrl_clean(pid_hcr_cor,pid_hsrl);
    
    % Combined by using both data sets in one process
    pid_comb2_cor=calc_pid_direct_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
        dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    
    %% Scales and units
    cscale_hsrl=[1,1,1;0,0,1.0;0,1,0;1,0.67,0;1,0,1;0,1,1;1,0.67,0];
    cscale_hcr=[1,1,1; 0,0,1.0; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0];
    cscale_comb=[1,1,1; 0,0,1; 0,1,0.; 1,0,0; 1,0,1; 0,1,1; 1,1,0; 0.5,0,0; 1,0.67,0];
    
    units_str_hsrl={'No signal','Cloud liquid','Drizzle',...
        'Aerosol1','SLW','Ice crystals','Aerosol2'};
    units_str_hcr={'No signal','Cloud liquid','Drizzle',...
        'Rain','SLW','Ice crystals','Snow','Wet snow/rimed ice'};
    units_str_comb={'No signal','Cloud liquid','Drizzle','Rain',...
        'SLW','Ice crystals','Snow','Wet snow/rimed ice','Aerosols'};
    
    %% Plot lidar
    close all
    if plotlidars==1
        plot_hsrl_clean(data,pid_hsrl,backscatLog,cscale_hsrl,units_str_hsrl,ylimits);
    end
    
    % Plot radar
    if plotradars==1
        plot_hcr_clean(data,pid_hcr_cor,cscale_hcr,units_str_hcr,ylimits);
    end
    
    % Plot radar and lidar
    if plotradars==1 & plotlidars==1
        plot_hsrl_hcr_clean(data,pid_comb_cor,backscatLog,cscale_comb,units_str_comb,ylimits);
    end
    
    %% PIDs
    % plot_pids_clean(data,pid_comb,pid_comb2,cscale_comb,units_str_comb,ylimits);
    
    plot_pids_clean(data,pid_comb_cor,pid_comb2_cor,cscale_comb,units_str_comb,ylimits);
        
end