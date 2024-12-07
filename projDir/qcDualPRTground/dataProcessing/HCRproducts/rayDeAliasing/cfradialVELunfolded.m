% add VEL_UNFOLDED data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz_combined'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';
whichModel='hrrr';
formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

%% Run processing

% Go through flights
for ii=5:size(caseList,1)
    
    disp(['IOP ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model=[];
        model.velUnfolded_long=[];
        
        model=read_model_longShort(model,modeldir,startTime,endTime);
        timeModelNum=datenum(model.time);
        
        % Check if times match
        if size(model.time,2)~=size(model.velUnfolded_long,2)
            error('Size of model time and model variable do not match.');
        end
        
        %% Loop through HCR data files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);

            % Check if variable exists
            try
                velIn=ncread(infile,'VEL_unfold_long');
            end
            if exist('velIn')
                warning('Variable already exists. Skipping file.')
                clear('velIn');
                continue
            end
            
            % Find times that are equal
            startTimeIn=ncread(infile,'time_coverage_start')';
            startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
                str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
            timeRead=ncread(infile,'time')';
            timeHCR=startTimeFile+seconds(timeRead);
            
            timeHcrNum=datenum(timeHCR);
            
            [C,ia,ib] = intersect(timeHcrNum,timeModelNum);
            
            if length(timeHCR)~=length(ib)
                warning('Times do not match up. Skipping file.')
                continue
            end
            
            % Write output
            fillVal=-9999;
            
            modVars=fields(model);
            
            for kk=1:length(modVars)
                if ~strcmp((modVars{kk}),'time')
                    modOut.(modVars{kk})=model.(modVars{kk})(:,ib);
                    modOut.(modVars{kk})(isnan(modOut.(modVars{kk})))=fillVal;
                    modOut.(modVars{kk})=modOut.(modVars{kk});
                end
            end
            
            % Open file
            ncid = netcdf.open(infile,'WRITE');
            netcdf.setFill(ncid,'FILL');
            
            % Get dimensions
            dimtime = netcdf.inqDimID(ncid,'time');
            dimrange = netcdf.inqDimID(ncid,'range');
            
            % Define variables
            netcdf.reDef(ncid);
            varidVEL = netcdf.defVar(ncid,'VEL_unfold_long','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidVEL,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidVEL,modOut.velUnfolded_long);
                       
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'VEL_unfold_long','long_name','doppler_velocity_unfolded');
            ncwriteatt(infile,'VEL_unfold_long','standard_name','radial_velocity_of_scatterers_away_from_instrument');
            ncwriteatt(infile,'VEL_unfold_long','units','m/s');
            ncwriteatt(infile,'VEL_unfold_long','comment','This field is created by de-aliasing velocity and masking non-cloud data.');
            ncwriteatt(infile,'VEL_unfold_long','grid_mapping','grid_mapping');
            ncwriteatt(infile,'VEL_unfold_long','coordinates','time range');
                        
        end
    end
end