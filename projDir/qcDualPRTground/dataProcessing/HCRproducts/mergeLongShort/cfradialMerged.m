% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz_combined'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';
whichModel='era5';

indir=HCRdir(project,quality,qcVersion,freqData);
[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    disp(['IOP ',num2str(ii)]);
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model=[];
        model.dbz=[];
        model.vel=[];
        model.width=[];
        model.ldr=[];
        
        model=read_model_longShort(model,modeldir,startTime,endTime);
        timeModelNum=datenum(model.time);
        
        %% Loop through HCR data files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);

            % Check if variable exists
            try
                dbzIn=ncread(infile,'DBZ');
            end
            if exist('dbzIn')
                warning('Variable already exists. Skipping file.')
                clear('dbzIn');
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
            fillVal=-999;
            
            modVars=fields(model);
            
            for kk=1:length(modVars)
                if ~strcmp((modVars{kk}),'time')
                    modOut.(modVars{kk})=model.(modVars{kk})(:,ib);
                    modOut.(modVars{kk})(isnan(modOut.(modVars{kk})))=fillVal;
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
            varidDBZ = netcdf.defVar(ncid,'DBZ','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidDBZ,false,fillVal);
            varidVEL = netcdf.defVar(ncid,'VEL','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidVEL,false,fillVal);
            varidWIDTH = netcdf.defVar(ncid,'WIDTH','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidWIDTH,false,fillVal);
            varidLDR = netcdf.defVar(ncid,'LDR','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidLDR,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidDBZ,modOut.dbz);
            netcdf.putVar(ncid,varidVEL,modOut.vel);
            netcdf.putVar(ncid,varidWIDTH,modOut.width);
            netcdf.putVar(ncid,varidLDR,modOut.ldr);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'DBZ','long_name','reflectivity');
            ncwriteatt(infile,'DBZ','standard_name','equivalent_reflectivity_factor');
            ncwriteatt(infile,'DBZ','units','dBZ');
            ncwriteatt(infile,'DBZ','comment','This field is created by merging the long and short pulse data.');
            ncwriteatt(infile,'DBZ','grid_mapping','grid_mapping');
            ncwriteatt(infile,'DBZ','coordinates','time range');

            ncwriteatt(infile,'VEL','long_name','doppler_velocity');
            ncwriteatt(infile,'VEL','standard_name','radial_velocity_of_scatterers_away_from_instrument');
            ncwriteatt(infile,'VEL','units','m/s');
            ncwriteatt(infile,'VEL','comment','This field is created by merging the long and short pulse data.');
            ncwriteatt(infile,'VEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'VEL','coordinates','time range');
                        
            ncwriteatt(infile,'WIDTH','long_name','spectrum_width');
            ncwriteatt(infile,'WIDTH','standard_name','doppler_spectrum_width');
            ncwriteatt(infile,'WIDTH','units','m/s');
            ncwriteatt(infile,'WIDTH','comment','This field is created by merging the long and short pulse data.');
            ncwriteatt(infile,'WIDTH','grid_mapping','grid_mapping');
            ncwriteatt(infile,'WIDTH','coordinates','time range');

            ncwriteatt(infile,'LDR','long_name','linear_depolarization_ratio');
            ncwriteatt(infile,'LDR','standard_name','log_linear_depolarization_ratio');
            ncwriteatt(infile,'LDR','units','dB');
            ncwriteatt(infile,'LDR','comment','This field is created by merging the long and short pulse data.');
            ncwriteatt(infile,'LDR','grid_mapping','grid_mapping');
            ncwriteatt(infile,'LDR','coordinates','time range');            
        end
    end
end