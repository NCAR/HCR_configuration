% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.1';
freqData='10hz'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

saveTime=0;

plotIn.plotMR=0;
plotIn.plotMax=0;

convThresh=4;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

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
    data.LDR=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    %data.SNR=[];

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
    data.LDR(isnan(data.DBZ_MASKED))=nan;
    
    ylimits=[0 (max(data.asl(~isnan(data.DBZ_MASKED)))./1000)+0.5];
    plotIn.ylimits=ylimits;

    tempOrig=data.TEMP;

    %% Calculate velocity texture

    pixRadVEL=50;
    velBase=-20;

    data.VEL_MASKED(:,data.elevation<0)=-data.VEL_MASKED(:,data.elevation<0);
    data.VELTEXT=f_velTexture(data.VEL_MASKED,data.elevation,pixRadVEL,velBase);
    
    %% Pre process

    disp('Pre processing ...');
    data=preProcessPID(data,convThresh);

    %% Calculate PID

    % HCR
    disp('Getting PID ...');
    pid_hcr=calc_pid(data,plotIn);

    %% Set areas above melting layer with no WIDTH and no LDR to cloud or precip

    data.TEMP=tempOrig;

    %     smallInds=find((data.MELTING_LAYER==20 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & (pid_hcr==3 | pid_hcr==6));
    %     pid_hcr(smallInds)=11;
    %
    %     largeInds=find((data.MELTING_LAYER==20 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
    %         (pid_hcr==1 | pid_hcr==2 | pid_hcr==4 | pid_hcr==5));
    %     pid_hcr(largeInds)=10;

    smallInds=find((data.MELTING_LAYER==20 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
        (data.DBZ_MASKED<=5 | data.VEL_MASKED<=1));
    pid_hcr(smallInds)=11;

    largeInds=find((data.MELTING_LAYER==20 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
        (data.DBZ_MASKED>5 & data.VEL_MASKED>1));
    pid_hcr(largeInds)=10;

    %% Set low DBZ to cloud liquid

    pid_hcr(pid_hcr>9 & data.DBZ_MASKED<=-30 & data.TEMP>-40)=3;

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