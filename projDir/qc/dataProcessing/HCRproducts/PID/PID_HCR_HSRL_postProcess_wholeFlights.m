% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, otrec, cset
whichModel='era5';

whichFilter=1; % 0: no filter, 1: mode filter, 2: coherence filter

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/combined/'];

%[~,directories.modeldir]=modelDir(project,whichModel,freqData);
%outdir=directories.modeldir;
outdir='/run/media/romatsch/RSF0006/rsf/pid_hcr/socratesMat/';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=3:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading data ...')
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
    
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
    
    disp('Creating HSRL PID ...');
    
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
        disp('Filtering ...');
        pid_hsrl=modeFilter(pid_hsrl,7,0.7);
    elseif whichFilter==2
        disp('Filtering ...');
        pid_hsrl=coherenceFilter(pid_hsrl,7,0.7);
    end
    
    %% Calculate HCR without attenuation correction
    
    disp('Creating HCR uncorrected PID ...');
    
    [pid_hcr]=calc_pid_hcr_postProcess(data.HCR_DBZ,data);
    pid_hcr(isnan(data.HCR_DBZ))=nan;
    
    % Combined from merging hcr and hsrl pid
    pid_comb=combine_pid_hcr_hsrl_postProcess(pid_hcr,pid_hsrl);
    
    %         % Combined by using both data sets in one process
    %         pid_comb2=calc_pid_direct_clean_eff(data.HSRL_Aerosol_Backscatter_Coefficient,lin_depol,...
    %             data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    
    %% Calculate attenuation correction
    
    disp('Calculating attenuation correction ...');
    
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
    
    disp('Creating corrected HCR PID ...');
    
    % HCR
    [pid_hcr_cor]=calc_pid_hcr_postProcess(dBZ_cor,data);
    pid_hcr_cor(isnan(dBZ_cor))=nan;
    
    if whichFilter==1
        disp('Filtering ...');
        pid_hcr_cor=modeFilter(pid_hcr_cor,7,0.7);
    elseif whichFilter==2
        disp('Filtering ...');
        pid_hcr_cor=coherenceFilter(pid_hcr_cor,7,0.7);
    end
    
    disp('Combining PIDs ...');
    
    % Combined from merging hcr and hsrl pid
    [pid_comb_cor ~]=combine_pid_hcr_hsrl_postProcess(pid_hcr_cor,pid_hsrl);
    
    %% Save
    disp('Saving PID fields ...')
    
    pidComb=pid_comb_cor;
    save([outdir,whichModel,'.pidComb.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pidComb');
    
    pidHCR=pid_hcr_cor;
    save([outdir,whichModel,'.pidHCR.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pidHCR');
    
    pidHSRL=pid_hsrl;
    save([outdir,whichModel,'.pidHSRL.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pidHSRL');
    
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
end