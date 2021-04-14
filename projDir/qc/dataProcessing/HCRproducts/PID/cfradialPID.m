% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; %socrates, aristo, cset

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%indir=HCRdir(project,quality,freqData);
indir='/run/media/romatsch/RSF0006/rsf/pid_hcr/socrates/';

%[~,modeldir]=modelDir(project,whichModel,freqData);
modeldir='/run/media/romatsch/RSF0006/rsf/pid_hcr/socratesMat/';

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model=[];
        model.pidComb=[];
        model.pidHCR=[];
        model.pidHSRL=[];
        
        model=read_model(model,modeldir,startTime,endTime);
        timeModelNum=datenum(model.time);
        
        %% Loop through HCR data files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);
            
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
            fillVal=-99;
            
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
            varidPID_COMB = netcdf.defVar(ncid,'PID_COMBINED','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidPID_COMB,false,fillVal);
            netcdf.endDef(ncid);
            
            netcdf.reDef(ncid);
            varidPID_HCR = netcdf.defVar(ncid,'PID_HCR','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidPID_HCR,false,fillVal);
            netcdf.endDef(ncid);
            
            netcdf.reDef(ncid);
            varidPID_HSRL = netcdf.defVar(ncid,'PID_HSRL','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidPID_HSRL,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidPID_COMB,modOut.pidComb);
            netcdf.putVar(ncid,varidPID_HCR,modOut.pidHCR);
            netcdf.putVar(ncid,varidPID_HSRL,modOut.pidHSRL);
                       
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'PID_COMBINED','long_name','particle_id_hcr_hsrl_combined');
            ncwriteatt(infile,'PID_COMBINED','standard_name','hydrometeor_type');
            ncwriteatt(infile,'PID_COMBINED','units','');
            ncwriteatt(infile,'PID_COMBINED','flag_values',[1, 2, 3, 4, 5, 6, 7, 8]);
            ncwriteatt(infile,'PID_COMBINED','flag_meanings','cloud_liquid drizzle rain slw ice_crystals snow wet_snow_rimed_ice aerosols');
            ncwriteatt(infile,'PID_COMBINED','grid_mapping','grid_mapping');
            ncwriteatt(infile,'PID_COMBINED','coordinates','time range');
            
            ncwriteatt(infile,'PID_HCR','long_name','particle_id_hcr');
            ncwriteatt(infile,'PID_HCR','standard_name','hydrometeor_type');
            ncwriteatt(infile,'PID_HCR','units','');
            ncwriteatt(infile,'PID_HCR','flag_values',[1, 2, 3, 4, 5, 6, 7]);
            ncwriteatt(infile,'PID_HCR','flag_meanings','cloud_liquid drizzle rain slw ice_crystals snow wet_snow_rimed_ice');
            ncwriteatt(infile,'PID_HCR','grid_mapping','grid_mapping');
            ncwriteatt(infile,'PID_HCR','coordinates','time range');
            
            ncwriteatt(infile,'PID_HSRL','long_name','particle_id_hsrl');
            ncwriteatt(infile,'PID_HSRL','standard_name','hydrometeor_type');
            ncwriteatt(infile,'PID_HSRL','units','');
            ncwriteatt(infile,'PID_HSRL','flag_values',[1, 2, 3, 4, 5, 6]);
            ncwriteatt(infile,'PID_HSRL','flag_meanings','cloud_liquid drizzle aerosols1 slw ice_crystals aerosols2');
            ncwriteatt(infile,'PID_HSRL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'PID_HSRL','coordinates','time range');
                        
        end
    end
end