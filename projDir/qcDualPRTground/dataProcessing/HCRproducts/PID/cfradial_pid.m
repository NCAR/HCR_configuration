% Add pid data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
qcVersion='v1.0';
freqData='10hz_combined';
whichModel='hrrr';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

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
        model.pid=[];
        
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
            varidPID_HCR = netcdf.defVar(ncid,'PID','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidPID_HCR,false,fillVal);
            netcdf.endDef(ncid);
                        
            % Write variables
            netcdf.putVar(ncid,varidPID_HCR,modOut.pid);
                                   
            netcdf.close(ncid);
            
            % Write attributes
                        
            ncwriteatt(infile,'PID','long_name','particle_id_hcr');
            ncwriteatt(infile,'PID','standard_name','hydrometeor_type');
            ncwriteatt(infile,'PID','units','');
            ncwriteatt(infile,'PID','flag_values',[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
            ncwriteatt(infile,'PID','flag_meanings',...
                'rain supercooled_rain drizzle supercooled_drizzle cloud_liquid supercooled_cloud_liquid melting large_frozen small_frozen precipitation cloud');
            ncwriteatt(infile,'PID','is_discrete','true');
            ncwriteatt(infile,'PID','grid_mapping','grid_mapping');
            ncwriteatt(infile,'PID','coordinates','time range');                                    
        end
    end
end
