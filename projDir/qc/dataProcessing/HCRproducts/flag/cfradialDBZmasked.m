% add model data to cfradial files
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='noreaster'; % socrates, cset, aristo, otrec
quality='qc2'; % field, qc1, qc2
qcVersion='v2.0';
freqData='10hz';
%whichModel='era5';

formatOut = 'yyyymmdd';

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

%% Run processing

% Go through flights
for ii=1:size(caseList,1)
    
    disp(['Flight ',num2str(ii)]);
    
    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        %% Loop through HCR data files
        for jj=1:length(fileList)
            
            testField=[];;
            
            infile=fileList{jj};
            disp(infile);
            
            try
                testField=ncread(infile,'DBZ_MASKED');
            end
            if ~isempty(testField)
                warning('Field already exists. Skipping file.')
                continue
            end                
                        
            % Create masked DBZ field
            dbz=ncread(infile,'DBZ')';
            mask=ncread(infile,'FLAG')';
            
            dbzMasked=dbz;
            dbzMasked(mask~=1)=nan;
            dbzMasked=dbzMasked';
            
            % Write output
            fillVal=-9999;
            
            % Open file
            ncid = netcdf.open(infile,'WRITE');
            netcdf.setFill(ncid,'FILL');
            
            % Get dimensions
            dimtime = netcdf.inqDimID(ncid,'time');
            dimrange = netcdf.inqDimID(ncid,'range');
            
            % Define variables
            netcdf.reDef(ncid);
            varidDBZm = netcdf.defVar(ncid,'DBZ_MASKED','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidDBZm,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidDBZm,dbzMasked);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'DBZ_MASKED','long_name','reflectivity');
            ncwriteatt(infile,'DBZ_MASKED','standard_name','equivalent_reflectivity_factor');
            ncwriteatt(infile,'DBZ_MASKED','units','dBZ');
            ncwriteatt(infile,'DBZ_MASKED','comment','This field is computed by applying DBZ(FLAG>1)=NAN');
            ncwriteatt(infile,'DBZ_MASKED','grid_mapping','grid_mapping');
            ncwriteatt(infile,'DBZ_MASKED','coordinates','time range');
            
        end
    end
end