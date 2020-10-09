% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc2'; %field, qc1, or qc2
freqData='2hzMerged'; % 10hz, 100hz, or 2hz
whichModel='era5';

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/combined_hcr_hsrl/',project,'/'];

%[~,directories.modeldir]=modelDir(project,whichModel,freqData);
%outdir=directories.modeldir;
outdir='/run/media/romatsch/RSF0006/rsf/pid_hcr/socratesMat/';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=6:size(caseList,1)
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
    
    disp('Calculating HSRL PID ...');
    
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
    
    disp('Calculating HCR PID without attenuation correction ...');
    
    [pid_hcr]=calc_pid_hcr_clean_eff(data.HCR_DBZ,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    pid_hcr(isnan(data.HCR_DBZ))=nan;
    pid_hcr(isnan(pid_hcr))=1;
    
    % Combined from merging hcr and hsrl pid
    pid_comb=combine_pid_hcr_hsrl_clean(pid_hcr,pid_hsrl);
   
    %% Calculate attenuation correction
    
    disp('Calculating attenuation correction ...');
    
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
    
    disp('Calculating HCR PID with attenuation correction ...');
    
    % HCR
    [pid_hcr_cor]=calc_pid_hcr_clean_eff(dBZ_cor,data.HCR_LDR,data.HCR_VEL,data.HCR_WIDTH,data.temp);
    pid_hcr_cor(isnan(dBZ_cor))=nan;
    pid_hcr_cor(isnan(pid_hcr_cor))=1;
    
    % Combined from merging hcr and hsrl pid
    pid_comb_cor=combine_pid_hcr_hsrl_clean(pid_hcr_cor,pid_hsrl); 
    
   %% Save
    disp('Saving PID ...')
    
    pid=pid_comb_cor;
    pid(pid==1)=nan;
    pid=pid-1;
    save([outdir,whichModel,'.pid.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pid');
    
    timeHCR=data.time;
    save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
end