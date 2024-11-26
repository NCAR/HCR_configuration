% Calculate PID from HCR HSRL combined data
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';
whichModel='hrrr';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

saveTime=0;

indir=HCRdir(project,quality,qcVersion,freqData);

[~,outdir]=modelDir(project,whichModel,quality,qcVersion,freqData);

caseList = table2array(readtable(infile));


plotIn.plotMR=0;
plotIn.plotMax=0;

convThresh=4;

whichFilter=0; % 0: no filter, 1: mode filter, 2: coherence filter
postProcess=1; % 1 if post processing is desired

for aa=1:size(caseList,1)
    disp(['IOP ',num2str(aa)]);
    disp('Loading data ...')

    startTime=datetime(caseList(aa,1:6));
    endTime=datetime(caseList(aa,7:12));

    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    disp([datestr(startTime,'yyyy-mm-dd HH:MM'),' to ',datestr(endTime,'yyyy-mm-dd HH:MM')]);

    %% Load data

    data=[];

    %HCR data
    data.DBZ=[];
    data.VEL=[];
    data.TEMP=[];
    data.MELTING_LAYER=[];
    %data.SNR=[];

    % Check if LDR_MASKED is available
    try
        velTest=ncread(fileList{1},'LDR_MASKED');
        data.LDR_MASKED=[];
    catch
        data.LDR=[];
    end

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    if isfield(data,'LDR_MASKED')
        data.LDR=data.LDR_MASKED;
        data=rmfield(data,'LDR_MASKED');
    end

    % Mask with FLAG
    data.LDR(isnan(data.DBZ))=nan;
    
    ylimits=[0 (max(data.asl(~isnan(data.DBZ)))./1000)+0.5];
    plotIn.ylimits=ylimits;

    tempOrig=data.TEMP;

    %% Calculate velocity texture

    pixRadVEL=300;
    velBase=-20;

    data.VEL(:,data.elevation<0)=-data.VEL(:,data.elevation<0);
    data.VELTEXT=f_velTexture(data.VEL,pixRadVEL,velBase);
    
    %% Pre process

    disp('Pre processing ...');
    data=preProcessPID(data,convThresh);

    %% Calculate PID

    % HCR
    disp('Getting PID ...');
    pid_hcr=calc_pid(data,plotIn);

    %% Set areas above melting layer with no WIDTH and no LDR to cloud or precip

    data.TEMP=tempOrig;

    smallInds=find((data.MELTING_LAYER>15 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
        (data.DBZ<=5 | data.VEL<=1));
    pid_hcr(smallInds)=11;

    largeInds=find((data.MELTING_LAYER>15 | (isnan(data.MELTING_LAYER) & data.TEMP<0)) & isnan(data.LDR) & ...
        (data.DBZ>5 & data.VEL>1));
    pid_hcr(largeInds)=10;

    %% Set low DBZ to cloud liquid

    pid_hcr(pid_hcr>9 & data.DBZ<=-30 & data.TEMP>-40)=3;

    %% Set Melting to melting
    pid_hcr(data.MELTING_LAYER==11 | data.MELTING_LAYER==19)=4;

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
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'pid','-v7.3');

    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.IOP',num2str(aa),'.mat'],'timeHCR');
    end

end