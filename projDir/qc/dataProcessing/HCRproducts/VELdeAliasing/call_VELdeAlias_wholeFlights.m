% De-alias velocity

clear all
close all

addpath(genpath('~/git/HCR_configuration/projDir/qc/dataProcessing/'));

project='socrates'; % socrates, cset, aristo, otrec
quality='qc3'; % field, qc1, qc2
qcVersion='v3.0';
freqData='10hz'; % 10hz, 100hz, or 2hz
whichModel='era5';

[~,directories.modeldir]=modelDir(project,whichModel,quality,qcVersion,freqData);

outdir=directories.modeldir;

infile=['~/git/HCR_configuration/projDir/qc/dataProcessing/scriptsFiles/flights_',project,'_data.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

for kk=1:size(caseList,1)
    
    disp(['Flight ',num2str(kk)]);
    disp('Loading HCR data.')
    
    startTime=datetime(caseList(kk,1:6));
    endTime=datetime(caseList(kk,7:12));
    
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
    
    if ~isempty(fileList)
        data=[];
        data.nyquist_velocity=[];
        data.VEL_CORR=[];
        data.FLAG=[];
               
        % Make list of files within the specified time frame
        fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
        
        if length(fileList)==0
            disp('No data files found.');
            return
        end
        
        % Load data
        data=read_HCR(fileList,data,startTime,endTime);
        
        velMasked=data.VEL_CORR;
        velMasked(data.FLAG~=1)=nan;
        
        %% Correct velocity folding
        
        disp('De-aliasing ...');
        velDeAliased=dealiasArea(velMasked,data.elevation,data.nyquist_velocity);
        
        %% Save
        disp('Saving velFinal field.')
        
        velUnfolded=velDeAliased;
        
        save([outdir,whichModel,'.velUnfolded.',datestr(data.time(1),'YYYYmmDD_HHMMSS'),'_to_',...
            datestr(data.time(end),'YYYYmmDD_HHMMSS'),'.Flight',num2str(kk),'.mat'],'velUnfolded');
        
    end
end