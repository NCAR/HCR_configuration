% add convstrat data to cfradial files

% Plot HCR convStrat from mat file in hourly plots

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%

project='socrates'; %socrates, aristo, cset, otrec
quality='qc3'; %field, qc1, or qc2
freqData='10hz';
qcVersion='v3.0';
whichModel='era5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%% Run processing

% Go through flights
for ii=2:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model.convStrat=[];
        model.convStrat1D=[];
        model.convectivity=[];
        
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
            varidConv = netcdf.defVar(ncid,'CONVECTIVITY','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidConv,false,fillVal);
            varidSC2D = netcdf.defVar(ncid,'ECHO_TYPE_2D','NC_SHORT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidSC2D,false,fillVal);
            varidSC1D = netcdf.defVar(ncid,'ECHO_TYPE_1D','NC_SHORT',[dimtime]);
            netcdf.defVarFill(ncid,varidSC1D,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidConv,modOut.convectivity);
            netcdf.putVar(ncid,varidSC2D,modOut.convStrat);
            netcdf.putVar(ncid,varidSC1D,modOut.convStrat1D);
                       
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'CONVECTIVITY','long_name','convective_probability');
            ncwriteatt(infile,'CONVECTIVITY','standard_name','convectivity');
            ncwriteatt(infile,'CONVECTIVITY','units','');
            ncwriteatt(infile,'CONVECTIVITY','grid_mapping','grid_mapping');
            ncwriteatt(infile,'CONVECTIVITY','coordinates','time range');
            
            ncwriteatt(infile,'ECHO_TYPE_2D','long_name','echo_type_2D');
            ncwriteatt(infile,'ECHO_TYPE_2D','standard_name','echo_type');
            ncwriteatt(infile,'ECHO_TYPE_2D','units','');
            ncwriteatt(infile,'ECHO_TYPE_2D','flag_values',[14, 16, 18, 25, 30, 32, 34, 36, 38]);
            ncwriteatt(infile,'ECHO_TYPE_2D','flag_meanings',...
                'strat_low strat_mid strat_high mixed conv conv_elevated conv_shallow conv_mid conv_deep');
            ncwriteatt(infile,'ECHO_TYPE_2D','is_discrete','true');
            ncwriteatt(infile,'ECHO_TYPE_2D','grid_mapping','grid_mapping');
            ncwriteatt(infile,'ECHO_TYPE_2D','coordinates','time range');
                        
            ncwriteatt(infile,'ECHO_TYPE_1D','long_name','echo_type_1D');
            ncwriteatt(infile,'ECHO_TYPE_1D','standard_name','echo_type');
            ncwriteatt(infile,'ECHO_TYPE_1D','units','');
            ncwriteatt(infile,'ECHO_TYPE_1D','flag_values',[14, 16, 18, 25, 30, 32, 34, 36, 38]);
            ncwriteatt(infile,'ECHO_TYPE_1D','flag_meanings',...
                'strat_low strat_mid strat_high mixed conv conv_elevated conv_shallow conv_mid conv_deep');
            ncwriteatt(infile,'ECHO_TYPE_1D','is_discrete','true');
            ncwriteatt(infile,'ECHO_TYPE_1D','grid_mapping','grid_mapping');
            ncwriteatt(infile,'ECHO_TYPE_1D','coordinates','time');
        end
    end
end
