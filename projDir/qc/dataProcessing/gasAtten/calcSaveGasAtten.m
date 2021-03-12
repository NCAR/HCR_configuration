% Calculate gaseous attenuation and save in cfradial files

% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='cset'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
freqData='10hz';
whichModel='era5';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

%indir=HCRdir(project,quality,freqData);
indir=['/run/media/romatsch/RSF0006/rsf/gasAtt/',project,'/10hz/'];

%% Run processing

% Go through flights
for ii=9:size(caseList,1)
    
    model=[];
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        
        
        %% Write cfradial files
        for jj=1:length(fileList)
            infile=fileList{jj};
            
            disp(infile);
            
            %             % Find times that are equal
            %             startTimeIn=ncread(infile,'time_coverage_start')';
            %             startTimeFile=datetime(str2num(startTimeIn(1:4)),str2num(startTimeIn(6:7)),str2num(startTimeIn(9:10)),...
            %                 str2num(startTimeIn(12:13)),str2num(startTimeIn(15:16)),str2num(startTimeIn(18:19)));
            %             timeRead=ncread(infile,'time')';
            %             timeHCR=startTimeFile+seconds(timeRead);
            %
            %             timeHcrNum=datenum(timeHCR);
            %
            %             [C,ia,ib] = intersect(timeHcrNum,timeModelNum);
            %
            %             if length(timeHCR)~=length(ib)
            %                 warning('Times do not match up. Skipping file.')
            %                 continue
            %             end
            
            data=[];
            
            data.DBZ=[];
            data.U_SURF=[];
            data.V_SURF=[];
            data.SST=[];
            data.TEMP=[];
            data.PRESS=[];
            data.RH=[];
            data.TOPO=[];
            data.FLAG=[];
            
            dataVars=fieldnames(data);
            
            % Load data
            data=read_HCR(fileList(jj),data,startTime,endTime);
            
            % Check if all variables were found
            for ii=1:length(dataVars)
                if ~isfield(data,dataVars{ii})
                    dataVars{ii}=[];
                end
            end
            
            data.frq=ncread(infile,'frequency');
            
            %% One way and two way gaseous attenuation
            
            disp('Calculating gaseous attenuation ...');
            [gasAttClear,gasAttCloud,gasAttClearMat,gasAttCloudMat]=get_gas_atten(data);
            gasAttCloud2=2*gasAttCloud';
            
            %% Prepare for writing
            model=[];
            model.specAtt=gasAttCloudMat;
            model.path2way=gasAttCloud2;
            
            model.time=data.time;
            timeModelNum=datenum(model.time);
            
            
            % Write output
            fillVal=-9999;
            
            modVars=fields(model);
            
            for kk=1:length(modVars)
                if ~strcmp((modVars{kk}),'time')
                    modOut.(modVars{kk})=model.(modVars{kk});
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
            varidSpec = netcdf.defVar(ncid,'SPECIFIC_GASEOUS_ATTENUATION','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidSpec,false,fillVal);
            varidPath = netcdf.defVar(ncid,'PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY','NC_FLOAT',[dimtime]);
            netcdf.defVarFill(ncid,varidPath,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidSpec,modOut.specAtt);
            netcdf.putVar(ncid,varidPath,modOut.path2way);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'SPECIFIC_GASEOUS_ATTENUATION','long_name','specific_gaseous_attenuation_1way');
            ncwriteatt(infile,'SPECIFIC_GASEOUS_ATTENUATION','standard_name','specific_gaseous_attenuation_1way');
            ncwriteatt(infile,'SPECIFIC_GASEOUS_ATTENUATION','units','dB');
            ncwriteatt(infile,'SPECIFIC_GASEOUS_ATTENUATION','grid_mapping','grid_mapping');
            ncwriteatt(infile,'SPECIFIC_GASEOUS_ATTENUATION','coordinates','time range');
            
            ncwriteatt(infile,'PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY','long_name','path_integrated_gaseous_attenuation_2way');
            ncwriteatt(infile,'PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY','standard_name','path_integrated_gaseous_attenuation_2way');
            ncwriteatt(infile,'PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY','units','dB');
            ncwriteatt(infile,'PATH_INTEGRATED_GASEOUS_ATTENUATION_2WAY','coordinates','time');
        end
    end
end