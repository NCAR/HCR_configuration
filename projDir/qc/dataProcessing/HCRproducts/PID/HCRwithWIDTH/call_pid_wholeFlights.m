% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

saveTime=0;

plotIn.plotMR=0;
plotIn.plotMax=0;

convThresh=4;
widthThresh=0.4;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

% if strcmp(project,'otrec')
%     indir='/scr/sleet2/rsfdata/projects/otrec/hcr/qc2/cfradial/development/pid/10hz/';
% elseif strcmp(project,'socrates')
%     indir='/scr/snow2/rsfdata/projects/socrates/hcr/qc2/cfradial/development/pid/10hz/';
% end
%
% outdir=[indir(1:end-30),'mat/pid/10hz/'];

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
    data.DBZ_MASKED=[];
    data.VEL_MASKED=[];
    data.WIDTH=[];
    data.LDR=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    data.SNR=[];

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
    data.WIDTH(isnan(data.DBZ_MASKED))=nan;
    data.LDR(isnan(data.DBZ_MASKED))=nan;
    data.TEMP(isnan(data.DBZ_MASKED))=nan;
    data.MELTING_LAYER(isnan(data.DBZ_MASKED))=nan;

    ylimits=[0 (max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5];
    plotIn.ylimits=ylimits;

    %% Correct for attenuation

    Z_lin=10.^(data.DBZ_MASKED*0.1);

    % Mask out non liquid data
    liqMeltInds=find(data.MELTING_LAYER<20);
    Z_lin(~liqMeltInds)=nan;

    wt_coef=nan(size(data.DBZ_MASKED));
    wt_exp=nan(size(data.DBZ_MASKED));

    wt_coef(data.DBZ_MASKED < - 20)=20.;
    wt_exp(data.DBZ_MASKED < - 20)=0.52;
    wt_coef(-20 <data.DBZ_MASKED <-15 )=1.73;
    wt_exp(-20 <data.DBZ_MASKED < -15 )=0.15;
    wt_coef(data.DBZ_MASKED > -15)=0.22;
    wt_exp(data.DBZ_MASKED > -15)=0.68;

    att_cumul=2.*0.0192*cumsum((wt_coef.*Z_lin.^wt_exp),1,'omitnan');
    dBZ_cor_all=data.DBZ_MASKED+att_cumul;

    % Replace dBZ values with attenuation corrected values in liquid and
    % melting regions
    dBZ_cor=data.DBZ_MASKED;
    dBZ_cor(liqMeltInds)=dBZ_cor_all(liqMeltInds);

    %% Calculate velocity texture

    pixRadVEL=50;
    velBase=-20;

    data.VELTEXT=f_velTexture(data.VEL_MASKED,data.elevation,pixRadVEL,velBase);

    %% Pre process

    disp('Pre processing ...');
    data=preProcessPID(data,convThresh,widthThresh);

    %% Calculate PID

    % HCR
    disp('Getting PID ...');

    %plotIn.figdir=[figdir,'debugPlots/'];

    pid_hcr=calc_pid(dBZ_cor,data,plotIn);

    %% Set areas above melting layer with no WIDTH and no LDR to cloud or precip

    smallInds=find(data.MELTING_LAYER==20 & isnan(data.LDR) & isnan(data.WIDTH) & (pid_hcr==3 | pid_hcr==6));
    pid_hcr(smallInds)=11;

    largeInds=find(data.MELTING_LAYER==20 & isnan(data.LDR) & isnan(data.WIDTH) & ...
        (pid_hcr==1 | pid_hcr==2 | pid_hcr==4 | pid_hcr==5));
    pid_hcr(largeInds)=10;

    %% Add supercooled

    disp('Adding supercooled ...')
    pid_hcr=addSupercooled(pid_hcr,data);

    %% Post process

    if postProcess
        disp('Post processing ...');
        pid_hcr=postProcessPID(pid_hcr,data);
    end

    %% Filter

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

    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end

end