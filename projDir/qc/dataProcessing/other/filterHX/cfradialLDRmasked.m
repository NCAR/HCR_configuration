% Created masked LDR field
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='spicule'; % socrates, cset, aristo, otrec
quality='qc1'; % field, qc1, qc2
qcVersion='v1.0';
freqData='10hz';
%whichModel='narr';

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
                testField=ncread(infile,'LDR_MASKED');
            end
            if ~isempty(testField)
                warning('Field already exists. Skipping file.')
                continue
            end
            
            
            
            data.DBMHX=ncread(infile,'DBMHX');
            data.LDR=ncread(infile,'LDR');
            data.FLAG=ncread(infile,'FLAG');
            
            data.LDR(data.FLAG>1)=nan;
            
            %% Find indices where correction needs to be applied
            medianHXend=median(data.DBMHX(760:770,:));
            
            corrInd=medianHXend>-96;
            
            if ~isempty(corrInd)
                dbm=data.DBMHX(:,corrInd==1);
                ldr=data.LDR(:,corrInd==1);
                
                %% Find filter indices
                med=medianHXend(corrInd==1);
                
                %dbmMask=dbm>-94;
                dbmMask=dbm>med+0.5;
                
                dbmMaskOut=zeros(size(dbmMask));
                
                for ii=1:size(dbmMask,2)
                    rayIn=dbmMask(:,ii);
                    dbmMaskOut(:,ii)=bwareaopen(rayIn,10);
                end
                
                maskUse=dbmMaskOut==0 & dbmMask==1;
                
                filterLDR=ldr;
                filterLDR(maskUse==1)=nan;
                
                LDR_masked=data.LDR;
                LDR_masked(:,corrInd==1)=filterLDR;
                
            else
                LDR_masked=data.LDR;
            end
            
            
            %
            %             % Create masked VEL field
            %             vel=ncread(infile,'VEL')';
            %             maskFlag=ncread(infile,'FLAG')';
            %             maskAnt=ncread(infile,'ANTFLAG')';
            %
            %             velMasked=vel;
            %             velMasked(maskAnt>2,:)=nan;
            %             velMasked(maskFlag>1)=nan;
            %             velMasked=velMasked';
            
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
            varidLDRm = netcdf.defVar(ncid,'LDR_MASKED','NC_FLOAT',[dimrange dimtime]);
            netcdf.defVarFill(ncid,varidLDRm,false,fillVal);
            netcdf.endDef(ncid);
            
            % Write variables
            netcdf.putVar(ncid,varidLDRm,LDR_masked);
            
            netcdf.close(ncid);
            
            % Write attributes
            ncwriteatt(infile,'LDR_MASKED','long_name','linear_depolarization_ratio');
            ncwriteatt(infile,'LDR_MASKED','standard_name','log_linear_depolarization_ratio');
            ncwriteatt(infile,'LDR_MASKED','units','dB');
            ncwriteatt(infile,'LDR_MASKED','comment','This field is created by masking out invalid and non-cloud data.');
            ncwriteatt(infile,'LDR_MASKED','grid_mapping','grid_mapping');
            ncwriteatt(infile,'LDR_MASKED','coordinates','time range');
            
        end
    end
end