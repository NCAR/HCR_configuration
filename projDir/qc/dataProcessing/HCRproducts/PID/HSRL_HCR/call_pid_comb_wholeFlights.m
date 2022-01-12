% Calculate PID from HCR HSRL combined data

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; %socrates, aristo, cset
quality='qc3'; %field, qc1, or qc2
qcVersion='v3.0';
freqData='combined'; % 10hz, 100hz, 2hz, or combined
whichModel='era5';

saveTime=1;

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
        
        disp('Loading data');
        
        data=[];
        
        % HCR
        data.HCR_DBZ=[];
        data.HCR_VEL=[];
        data.HCR_LDR=[];
        data.HCR_MELTING_LAYER=[];

        % HSRL
        data.HSRL_Aerosol_Backscatter_Coefficient=[];
        data.HSRL_Particle_Linear_Depolarization_Ratio=[];

        % TEMP
        data.TEMP=[];
                       
        dataVars=fieldnames(data);
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        % Check if all variables were found
        for ii=1:length(dataVars)
            if ~isfield(data,dataVars{ii})
                dataVars{ii}=[];
            end
        end
       
        ylimits=[0 (max(data.asl(~isnan(data.HCR_DBZ)))./1000)+0.5];
        plotIn.ylimits=ylimits;

        tempOrig=data.TEMP;

        %% Calculate velocity texture

        pixRadVEL=10;
        velBase=-20;

        data.VELTEXT=f_velTexture(data.HCR_VEL,data.elevation,pixRadVEL,velBase);

        %% Mask LDR and HSRL

        data.HCR_LDR(isnan(data.HCR_DBZ))=nan;
        data.HSRL_Aerosol_Backscatter_Coefficient(isnan(data.HCR_DBZ))=nan;
        data.HSRL_Particle_Linear_Depolarization_Ratio(isnan(data.HCR_DBZ))=nan;

        %% Pre process

        disp('Pre processing ...');
        data=preProcessPIDcomb(data,convThresh);

        %% Calculate PID

        % HCR
        disp('Getting PID ...');

        plotIn.figdir=[];

        pid_hcr_hsrl=calc_pid_comb(data,plotIn);

        %% Set areas above melting layer with no LDR to cloud or precip

        data.TEMP=tempOrig;

        smallInds=find((data.HCR_MELTING_LAYER==20 | isnan(data.HCR_MELTING_LAYER & data.TEMP<0)) & isnan(data.HCR_LDR) & isnan(data.HSRL_Particle_Linear_Depolarization_Ratio) & ...
            (pid_hcr_hsrl==3 | pid_hcr_hsrl==6));
        pid_hcr_hsrl(smallInds)=11;

        largeInds=find(data.HCR_MELTING_LAYER==20 | isnan(data.HCR_MELTING_LAYER & data.TEMP<0)) & isnan(data.HCR_LDR) & isnan(data.HSRL_Particle_Linear_Depolarization_Ratio) & ...
            (pid_hcr_hsrl==1 | pid_hcr_hsrl==2 | pid_hcr_hsrl==4 | pid_hcr_hsrl==5));
        pid_hcr_hsrl(largeInds)=10;

        %% Set low DBZ to cloud liquid

        pid_hcr_hsrl(pid_hcr_hsrl>9 & data.HCR_DBZ<=-30)=3;

        %% Add supercooled

        disp('Adding supercooled ...')
        pid_hcr_hsrl=addSupercooledComb(pid_hcr_hsrl,data);

        %% Post process

        if postProcess
            disp('Post processing ...');
            pid_hcr_hsrl=postProcessPIDcomb(pid_hcr_hsrl,data);
        end

        %% Filter

        if whichFilter==1
            pid_hcr_hsrl=modeFilter(pid_hcr_hsrl,7,0.7);
        elseif whichFilter==2
            pid_hcr_hsrl=coherenceFilter(pid_hcr_hsrl,7,0.7);
        end

    %% Save
    disp('Saving PID ...')

    pid=pid_hcr_hsrl;
    save([outdir,whichModel,'.pid.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
        datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'pid');

    if saveTime
        timeHCR=data.time;
        save([outdir,whichModel,'.time.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(aa),'.mat'],'timeHCR');
    end

end