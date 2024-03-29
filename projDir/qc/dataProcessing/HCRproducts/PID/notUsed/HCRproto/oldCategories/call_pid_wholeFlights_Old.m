% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; %socrates, aristo, cset
quality='qc0'; %field, qc1, or qc2
qcVersion='v0.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
whichModel='ecmwf';

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=0; % 1 if post processing is desired

indir=HCRdir(project,quality,qcVersion,freqData);
%indir=['/run/media/romatsch/RSF0006/rsf/meltingLayer/',project,'/combined/'];

[~,directories.modeldir]=modelDir(project,whichModel,freqData);
outdir=directories.modeldir;
%outdir='/run/media/romatsch/RSF0006/rsf/pid_hcr/socratesMat/';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

for aa=1:size(caseList,1)
    disp(['Flight ',num2str(aa)]);
    disp('Loading data ...')
    
    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);
    
    %% Load data
    
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
        disp('Creating PID ...');
        % HCR
        [pid_hcr]=calc_pid(dBZ_cor,data,postProcess);
        
        if whichFilter==1
            pid_hcr=modeFilter(pid_hcr,7,0.7);
        elseif whichFilter==2
            pid_hcr=coherenceFilter(pid_hcr,7,0.7);
        end
 
    %% Save
    disp('Saving PID ...')
    
    pid=pid_hcr;
    save([outdir,whichModel,'.pid.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pid');
 
end