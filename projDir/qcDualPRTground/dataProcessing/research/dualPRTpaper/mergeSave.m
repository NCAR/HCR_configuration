% find minimum reflectivity values
clear all;
close all;

addpath(genpath('~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/'));

project='meow';
quality='qc1';
freqData='10hz_combined';
qcVersion='v1.0';

infile=['~/git/HCR_configuration/projDir/qcDualPRTground/dataProcessing/scriptsFiles/iops_',project,'.txt'];

caseList = table2array(readtable(infile));

indir=HCRdir(project,quality,qcVersion,freqData);

figdir='/scr/virga1/rsfdata/projects/meow/hcr/qc1/cfradial/v1.0_full/dualPRTpaper/';

%% Run processing

% Go through iops
for ii=1:size(caseList,1)

    disp(['IOP ',num2str(ii),' of ',num2str(size(caseList,1))]);

    startTime=datetime(caseList(ii,1:6));
    endTime=datetime(caseList(ii,7:12));

    data=[];

    data.DBZ_short=[];
    data.VEL_unfold_short=[];
    data.WIDTH_short=[];
    data.SNRVC_short=[];
    data.LDRV_short=[];
    data.FLAG_short=[];
    data.DBZ_long=[];
    data.VEL_unfold_long=[];
    data.WIDTH_long=[];
    data.LDRV_long=[];
       
    %% Load data
    disp('Loading data ...');

    % Make list of files within the specified time frame
    fileList=makeFileList(indir,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

    % Load data
    data=read_HCR(fileList,data,startTime,endTime);

    data.WIDTH_long(data.WIDTH_long<=0.1 | data.FLAG_short~=1)=nan;
    data.WIDTH_short(data.WIDTH_short<=0.1 | data.FLAG_short~=1)=nan;
    data.DBZ_long(data.DBZ_long<=-40 | data.FLAG_short~=1)=nan;
    data.DBZ_short(data.DBZ_short<=-40 | data.FLAG_short~=1)=nan;
    data.LDRV_long(data.LDRV_long>-1 | data.FLAG_short~=1)=nan;
    data.LDRV_short(data.LDRV_short>-1 | data.FLAG_short~=1)=nan;
    data.VEL_long=data.VEL_unfold_long;
    data.VEL_short=data.VEL_unfold_short;
    data.VEL_long(data.FLAG_short~=1)=nan;
    data.VEL_short(data.FLAG_short~=1)=nan;

    %% Sort by SNR
    disp('Sorting by SNR ...');
    snrEdges.DBZ=-6:-5;
    snrEdges.VEL=-5:-4;
    snrEdges.WIDTH=10:11;
    snrEdges.LDRV=25:26;

    dataFields={'DBZ','VEL','WIDTH','LDRV'};

    bySNR=[];
    for jj=1:length(snrEdges.DBZ)-1
        for kk=1:size(dataFields,2)
            thisName=dataFields{kk};
            thisMask=(data.SNRVC_short>snrEdges.(thisName)(jj) & data.SNRVC_short<=snrEdges.(thisName)(jj+1));
            longField=data.([thisName,'_long']);
            shortField=data.([thisName,'_short']);
            shortLong=cat(2,shortField(thisMask==1),longField(thisMask==1));
            shortLong(any(isnan(shortLong),2),:)=[];
            bySNR.(['bin',num2str(jj)]).(thisName)=shortLong;
            if ii>1
                bySNRall.(['bin',num2str(jj)]).(thisName)=cat(1,bySNRall.(['bin',num2str(jj)]).(thisName),shortLong);
            end
        end
    end

    if ii==1
        bySNRall=bySNR;
    end

end

save([figdir,'mergeScatterCensor.mat'],'bySNRall');