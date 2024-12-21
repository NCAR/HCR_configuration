% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
qcVersion='v2.0';
freqData='10hz';
whichModel='era5';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

[~,modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

%% Run processing

% Go through flights
for ii=2:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        % Get model data
        model=[];
        model.momentsSpecParams=[];
                
        model=read_model(model,modeldir,startTime,endTime);
        model=rmfield(model,'momentsSpecParams');
        timeModelNum=datenum(model.time);
        
        %model.sst(model.topo>0)=nan;
        
        %% Loop through HCR data files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);

            % Check if variable exists
            try
                sk=ncread(infile,'SKEWNESS');
            end
            if exist('sk')
                warning('Variable already exists. Skipping file.')
                clear('sk');
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
                if ~strcmp((modVars{kk}),'time') & ~strcmp((modVars{kk}),'asl')
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
            varidWidth = netcdf.defVar(ncid,'WIDTH_SPEC','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidWidth,false,fillVal);
            varidSkew = netcdf.defVar(ncid,'SKEWNESS','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidSkew,false,fillVal);
            varidKurt = netcdf.defVar(ncid,'KURTOSIS','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidKurt,false,fillVal);
            varidLS = netcdf.defVar(ncid,'LEFT_SLOPE','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidLS,false,fillVal);
            varidRS = netcdf.defVar(ncid,'RIGHT_SLOPE','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidRS,false,fillVal);
            varidEEw = netcdf.defVar(ncid,'EDGE_EDGE_WIDTH','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidEEw,false,fillVal);
            varidLEv = netcdf.defVar(ncid,'LEFT_EDGE_VEL','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidLEv,false,fillVal);
            varidREv = netcdf.defVar(ncid,'RIGHT_EDGE_VEL','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidREv,false,fillVal);
            varidLPv = netcdf.defVar(ncid,'LEFT_PEAK_VEL','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidLPv,false,fillVal);
            varidRPv = netcdf.defVar(ncid,'RIGHT_PEAK_VEL','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidRPv,false,fillVal);
           
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidWidth,modOut.width);
            netcdf.putVar(ncid,varidSkew,modOut.skew);
            netcdf.putVar(ncid,varidKurt,modOut.kurt);
            netcdf.putVar(ncid,varidLS,modOut.lslope);
            netcdf.putVar(ncid,varidRS,modOut.rslope);
            netcdf.putVar(ncid,varidEEw,modOut.lrwidth);
            netcdf.putVar(ncid,varidLEv,modOut.level);
            netcdf.putVar(ncid,varidREv,modOut.revel);
            netcdf.putVar(ncid,varidLPv,modOut.lpvel);
            netcdf.putVar(ncid,varidRPv,modOut.rpvel);
                        
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'WIDTH_SPEC','long_name','spectrum_width_from_spectra');
            ncwriteatt(infile,'WIDTH_SPEC','standard_name','doppler_spectrum_width');
            ncwriteatt(infile,'WIDTH_SPEC','units','m/s');
            ncwriteatt(infile,'WIDTH_SPEC','grid_mapping','grid_mapping');
            ncwriteatt(infile,'WIDTH_SPEC','coordinates','time range');
            
            ncwriteatt(infile,'SKEWNESS','long_name','spectral_skewness');
            ncwriteatt(infile,'SKEWNESS','standard_name','doppler_spectrum_skewness');
            ncwriteatt(infile,'SKEWNESS','units','m/s');
            ncwriteatt(infile,'SKEWNESS','grid_mapping','grid_mapping');
            ncwriteatt(infile,'SKEWNESS','coordinates','time range');
            
            ncwriteatt(infile,'KURTOSIS','long_name','spectral_kurtosis');
            ncwriteatt(infile,'KURTOSIS','standard_name','doppler_spectrum_kurtosis');
            ncwriteatt(infile,'KURTOSIS','units','m/s');
            ncwriteatt(infile,'KURTOSIS','grid_mapping','grid_mapping');
            ncwriteatt(infile,'KURTOSIS','coordinates','time range');

            ncwriteatt(infile,'LEFT_SLOPE','long_name','left_slope_of_spectrum');
            ncwriteatt(infile,'LEFT_SLOPE','standard_name','doppler_spectrum_left_slope');
            ncwriteatt(infile,'LEFT_SLOPE','units','(dB s)/m');
            ncwriteatt(infile,'LEFT_SLOPE','grid_mapping','grid_mapping');
            ncwriteatt(infile,'LEFT_SLOPE','coordinates','time range');
            
            ncwriteatt(infile,'RIGHT_SLOPE','long_name','right_slope_of_spectrum');
            ncwriteatt(infile,'RIGHT_SLOPE','standard_name','doppler_spectrum_right_slope');
            ncwriteatt(infile,'RIGHT_SLOPE','units','(dB s)/m');
            ncwriteatt(infile,'RIGHT_SLOPE','grid_mapping','grid_mapping');
            ncwriteatt(infile,'RIGHT_SLOPE','coordinates','time range');

            ncwriteatt(infile,'EDGE_EDGE_WIDTH','long_name','edge_to_edge_width_of_spectrum');
            ncwriteatt(infile,'EDGE_EDGE_WIDTH','standard_name','doppler_spectrum_edge_to_edge_width');
            ncwriteatt(infile,'EDGE_EDGE_WIDTH','units','m/s');
            ncwriteatt(infile,'EDGE_EDGE_WIDTH','grid_mapping','grid_mapping');
            ncwriteatt(infile,'EDGE_EDGE_WIDTH','coordinates','time range');

            ncwriteatt(infile,'LEFT_EDGE_VEL','long_name','left_edge_velocity_of_spectrum');
            ncwriteatt(infile,'LEFT_EDGE_VEL','standard_name','doppler_spectrum_left_edge_velocity');
            ncwriteatt(infile,'LEFT_EDGE_VEL','units','m/s');
            ncwriteatt(infile,'LEFT_EDGE_VEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'LEFT_EDGE_VEL','coordinates','time range');

            ncwriteatt(infile,'RIGHT_EDGE_VEL','long_name','right_edge_velocity_of_spectrum');
            ncwriteatt(infile,'RIGHT_EDGE_VEL','standard_name','doppler_spectrum_right_edge_velocity');
            ncwriteatt(infile,'RIGHT_EDGE_VEL','units','m/s');
            ncwriteatt(infile,'RIGHT_EDGE_VEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'RIGHT_EDGE_VEL','coordinates','time range');

            ncwriteatt(infile,'LEFT_PEAK_VEL','long_name','left_peak_velocity_of_spectrum');
            ncwriteatt(infile,'LEFT_PEAK_VEL','standard_name','doppler_spectrum_left_peak_velocity');
            ncwriteatt(infile,'LEFT_PEAK_VEL','units','m/s');
            ncwriteatt(infile,'LEFT_PEAK_VEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'LEFT_PEAK_VEL','coordinates','time range');

            ncwriteatt(infile,'RIGHT_PEAK_VEL','long_name','right_peak_velocity_of_spectrum');
            ncwriteatt(infile,'RIGHT_PEAK_VEL','standard_name','doppler_spectrum_right_peak_velocity');
            ncwriteatt(infile,'RIGHT_PEAK_VEL','units','m/s');
            ncwriteatt(infile,'RIGHT_PEAK_VEL','grid_mapping','grid_mapping');
            ncwriteatt(infile,'RIGHT_PEAK_VEL','coordinates','time range');
            
        end
    end
end