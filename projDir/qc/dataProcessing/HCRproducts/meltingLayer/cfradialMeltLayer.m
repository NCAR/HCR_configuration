% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='otrec'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz';
whichModel='era5';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,freqData);

[~,modeldir]=modelDir(project,whichModel,freqData);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model.meltLayer=[];
        
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
            varidML = netcdf.defVar(ncid,'FREEZING_LEVEL','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidML,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidML,modOut.meltLayer);
                       
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'FREEZING_LEVEL','long_name','freezing_level_and_zero_degree_level');
            ncwriteatt(infile,'FREEZING_LEVEL','standard_name','freezing_level_and_zero_degree_level');
            ncwriteatt(infile,'FREEZING_LEVEL','units','');
            ncwriteatt(infile,'FREEZING_LEVEL','flag_values',[0, 1, 2, 3]);
            ncwriteatt(infile,'FREEZING_LEVEL','flag_meanings','ERA5_zero_degree_level freezing_level_detected freezing_level_interpolated freezing_level_estimated');
            ncwriteatt(infile,'FREEZING_LEVEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'FREEZING_LEVEL','coordinates','time range');
                        
        end
    end
end