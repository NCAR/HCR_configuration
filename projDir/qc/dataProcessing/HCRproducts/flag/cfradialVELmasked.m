% Created masked velocity field
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2
qcVersion='v3.2';
freqData='10hz';

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
                testField=ncread(infile,'VEL_MASKED');
            end
            if ~isempty(testField)
                warning('Field already exists. Skipping file.')
                continue
            end                
                        
            % Create masked VEL field
            vel=ncread(infile,'VEL_UNFOLDED')';
            % maskFlag=ncread(infile,'FLAG')';
            maskAnt=ncread(infile,'ANTFLAG')';
            
            velMasked=vel;
            velMasked(maskAnt>2,:)=nan;
            % velMasked(maskFlag~=1)=nan;
            velMasked=velMasked';
            
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
            varidVELm = netcdf.defVar(ncid,'VEL_MASKED','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidVELm,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidVELm,velMasked);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'VEL_MASKED','long_name','doppler_velocity');
            ncwriteatt(infile,'VEL_MASKED','standard_name','radial_velocity_of_scatterers_away_from_instrument');
            ncwriteatt(infile,'VEL_MASKED','units','m/s');
            ncwriteatt(infile,'VEL_MASKED','comment','This field is created by applying VEL_UNFOLDED(ANTFLAG>2)=NAN');
            ncwriteatt(infile,'VEL_MASKED','grid_mapping','grid_mapping');
            ncwriteatt(infile,'VEL_MASKED','coordinates','time range');
            
        end
    end
end