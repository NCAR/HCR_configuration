% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow'; %socrates, aristo, cset
quality='qc1'; %field, qc1, or qc2
freqData='10hz_combined'; % 10hz, 100hz, or 2hz
qcVersion='v1.0';
whichModel='hrrr';

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
        model.antstat_short=[];
        model.antstat_long=[];
        model.flagfield_short=[];
        model.flagfield_long=[];
        
        model=read_model_longShort(model,modeldir,startTime,endTime);
        if length(model.antstat_long)~=length(model.antstat_short) | ...
                size(model.flagfield_long,2)~=size(model.flagfield_short,2)
            error('Short and long data sizes do not agree.')
        end
        timeModelNum=datenum(model.time);
        
        %% Loop through HCR data files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);

            % Check if variable exists
            try
                flagIn=ncread(infile,'FLAG_long');
            end
            if exist('flagIn')
                warning('Variable already exists. Skipping file.')
                clear('flagIn');
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
            fillVal=-99;
            
            modVars=fields(model);
            
            for kk=1:length(modVars)
                if ~strcmp((modVars{kk}),'time')
                    modOut.(modVars{kk})=model.(modVars{kk})(:,ib);
                    modOut.(modVars{kk})(isnan(modOut.(modVars{kk})))=fillVal;
                    modOut.(modVars{kk})=int16(modOut.(modVars{kk}));
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
            varidFLAGs = netcdf.defVar(ncid,'FLAG_short','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidFLAGs,false,fillVal);
            varidFLAGl = netcdf.defVar(ncid,'FLAG_long','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidFLAGl,false,fillVal);
            varidANTs = netcdf.defVar(ncid,'ANTFLAG_short','NC_SHORT',[dimtime]);
            netcdf.defVarFill(ncid,varidANTs,false,fillVal);
            varidANTl = netcdf.defVar(ncid,'ANTFLAG_long','NC_SHORT',[dimtime]);
            netcdf.defVarFill(ncid,varidANTl,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidFLAGs,modOut.flagfield_short);
            netcdf.putVar(ncid,varidFLAGl,modOut.flagfield_long);
            netcdf.putVar(ncid,varidANTs,modOut.antstat_short);
            netcdf.putVar(ncid,varidANTl,modOut.antstat_long);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'FLAG_short','long_name','data_flag');
            ncwriteatt(infile,'FLAG_short','standard_name','data_flag');
            ncwriteatt(infile,'FLAG_short','units','');
            ncwriteatt(infile,'FLAG_short','flag_values',[1,2,3,4,5,6,7,8,9,10,11]);
            ncwriteatt(infile,'FLAG_short','flag_meanings',...
                'cloud speckle extinct backlobe out_of_range transmitter_pulse water_surface land_surface below_surface noise_source_cal missing');
            ncwriteatt(infile,'FLAG_short','is_discrete','true');
            ncwriteatt(infile,'FLAG_short','grid_mapping','grid_mapping');
            ncwriteatt(infile,'FLAG_short','coordinates','time range');

            ncwriteatt(infile,'FLAG_long','long_name','data_flag');
            ncwriteatt(infile,'FLAG_long','standard_name','data_flag');
            ncwriteatt(infile,'FLAG_long','units','');
            ncwriteatt(infile,'FLAG_long','flag_values',[1,2,3,4,5,6,7,8,9,10,11]);
            ncwriteatt(infile,'FLAG_long','flag_meanings',...
                'cloud speckle extinct backlobe out_of_range transmitter_pulse water_surface land_surface below_surface noise_source_cal missing');
            ncwriteatt(infile,'FLAG_long','is_discrete','true');
            ncwriteatt(infile,'FLAG_long','grid_mapping','grid_mapping');
            ncwriteatt(infile,'FLAG_long','coordinates','time range');
                        
            ncwriteatt(infile,'ANTFLAG_short','long_name','antenna_flag');
            ncwriteatt(infile,'ANTFLAG_short','standard_name','antenna_flag');
            ncwriteatt(infile,'ANTFLAG_short','units','');
            ncwriteatt(infile,'ANTFLAG_short','flag_values',[1,2,3,4,5,6]);
            ncwriteatt(infile,'ANTFLAG_short','flag_meanings',...
                'down up pointing scanning transition failure');
            ncwriteatt(infile,'ANTFLAG_short','is_discrete','true');
            ncwriteatt(infile,'ANTFLAG_short','coordinates','time');

            ncwriteatt(infile,'ANTFLAG_long','long_name','antenna_flag');
            ncwriteatt(infile,'ANTFLAG_long','standard_name','antenna_flag');
            ncwriteatt(infile,'ANTFLAG_long','units','');
            ncwriteatt(infile,'ANTFLAG_long','flag_values',[1,2,3,4,5,6]);
            ncwriteatt(infile,'ANTFLAG_long','flag_meanings',...
                'down up pointing scanning transition failure');
            ncwriteatt(infile,'ANTFLAG_long','is_discrete','true');
            ncwriteatt(infile,'ANTFLAG_long','coordinates','time');
            
        end
    end
end